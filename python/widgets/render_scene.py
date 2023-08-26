from PySide6.QtOpenGLWidgets import QOpenGLWidget
from PySide6.QtWidgets import QApplication, QMessageBox 
from PySide6.QtCore import QSize, QRect, QTimer, Qt, Signal, Slot, QPointF
from PySide6.QtOpenGL import QOpenGLVertexArrayObject, QOpenGLBuffer, \
                             QOpenGLShaderProgram, QOpenGLShader
from PySide6.QtGui import QVector3D, QOpenGLFunctions, QMatrix4x4, \
                            QOpenGLContext, QSurfaceFormat, QVector3DList, QWheelEvent
import numpy as np
from shiboken6 import VoidPtr
import ctypes
import math
import sys

sys.path.append("C:\\Users\\Student\\OneDrive - Bournemouth University\\Desktop\\Personal\\Python\\Fluid_Predictor\\python\\")

from Fluid_Utilities.system import FluidSystem

try:
    from OpenGL import GL

except ImportError:
    app = QApplication(sys.argv)
    message_box = QMessageBox(QMessageBox.Critical, "OpenGL hellogl",
                                         "PyOpenGL must be installed to run this example.",
                                         QMessageBox.Close)
    message_box.setDetailedText("Run:\npip install PyOpenGL PyOpenGL_accelerate")
    message_box.exec()
    sys.exit(1)

class FluidSolver:

    def __init__(self, system_obj:FluidSystem):
        
        self.system_obj = system_obj
        self.m_data = QVector3DList()
        self.m_data.reserve(system_obj.num_particles)
        self.initialize_m_data()

    def initialize_m_data(self):
        for particle in self.system_obj.particle_list:
            self.add(QVector3D(particle.initial_pos[0],
                               particle.initial_pos[1],
                               particle.initial_pos[2]))
            
    def update_m_data(self, frame_index):
        self.m_data.clear()
        for i in range(self.system_obj.num_particles):
            position = self.system_obj.stored_positions[frame_index][i]
            self.add(QVector3D(position[0], position[1], position[2])) 

    def const_data(self):
        return self.m_data.constData()

    def count(self):
        return len(self.m_data) * 3

    def vertex_count(self):
        return self.count() / 6

    def add(self, v):
        self.m_data.append(v)


class RenderScene(QOpenGLWidget, QOpenGLFunctions):
    x_rotation_changed = Signal(int)
    y_rotation_changed = Signal(int)
    z_rotation_changed = Signal(int)

    x_pan_changed = Signal(int)
    y_pan_changed = Signal(int)

    def __init__(self, transparent, 
                 system_obj:FluidSystem, parent=None):
        
        QOpenGLWidget.__init__(self, parent)
        QOpenGLFunctions.__init__(self)

        self._transparent = transparent
        self._core = QSurfaceFormat.defaultFormat().profile() == QSurfaceFormat.CoreProfile

        self.system_obj = system_obj

        self._x_rot = 0
        self._y_rot = 0
        self._z_rot = 0

        self._x_pan = 0
        self._y_pan = 0

        self.zoom = 1.0
        self._last_pos = QPointF()
        self.fluid_solver = FluidSolver(self.system_obj)

        self.current_frame_index = 0

        self.all_frames_updated = False
        self.animation_timer = QTimer()
        self.animation_timer.timeout.connect(self.animate)
        self.animation_timer.start(1000 // 30)

        self.vao = QOpenGLVertexArrayObject()
        self._solver_vbo = QOpenGLBuffer()
        self.program = QOpenGLShaderProgram()
        self.fbo_shader_program = QOpenGLShaderProgram()
        self._proj_matrix_loc = 0
        self._mv_matrix_loc = 0
        self._normal_matrix_loc = 0
        self._light_pos_loc = 0
        self.proj = QMatrix4x4()
        self.camera = QMatrix4x4()
        self.world = QMatrix4x4()
        if transparent:
            fmt = self.format()
            fmt.setAlphaBufferSize(8)
            self.setFormat(fmt)

        """ self.system_obj.update_stored_positions() """

    def x_rotation(self):
        return self._x_rot

    def y_rotation(self):
        return self._y_rot

    def z_rotation(self):
        return self._z_rot

    def x_pan(self):
        return self._x_pan

    def y_pan(self):
        return self._y_pan
    
    def minimumSizeHint(self):
        return QSize(50, 50)

    def sizeHint(self):
        return QSize(400, 400)

    def normalize_angle(self, angle):
        while angle < 0:
            angle += 360 * 16
        while angle > 360 * 16:
            angle -= 360 * 16
        return angle

    @Slot(int)
    def set_xrotation(self, angle):
        angle = self.normalize_angle(angle)
        if angle != self._x_rot:
            self._x_rot = angle
            self.x_rotation_changed.emit(angle)
            self.update()

    @Slot(int)
    def set_yrotation(self, angle):
        angle = self.normalize_angle(angle)
        if angle != self._y_rot:
            self._y_rot = angle
            self.y_rotation_changed.emit(angle)
            self.update()

    @Slot(int)
    def set_zrotation(self, angle):
        angle = self.normalize_angle(angle)
        if angle != self._z_rot:
            self._z_rot = angle
            self.z_rotation_changed.emit(angle)
            self.update()

    @Slot(int)
    def set_xpan(self, offset):
        if offset != self._x_pan:
            self._x_pan = offset
            self.x_pan_changed.emit(offset)
            self.update()

    @Slot(int)
    def set_ypan(self, offset):
        if offset != self._y_pan:
            self._y_pan = offset
            self.y_pan_changed.emit(offset)
            self.update()

    @Slot()
    def cleanup(self):
        self.makeCurrent()
        self._logo_vbo.destroy()
        del self.program
        self.program = None
        self.doneCurrent()

    def vertex_shader_source_core(self):
        return """#version 150
                in vec4 vertex;
                in vec3 normal;
                out vec3 vert;
                out vec3 vertNormal;
                uniform mat4 projMatrix;
                uniform mat4 mvMatrix;
                uniform mat3 normalMatrix;
                void main() {
                   vert = vertex.xyz;
                   vertNormal = normalMatrix * normal;
                   gl_Position = projMatrix * mvMatrix * vertex;
                }"""

    def fragment_shader_source_core(self):
        return """#version 150
                in highp vec3 vert;
                in highp vec3 vertNormal;
                out highp vec4 fragColor;
                uniform highp vec3 lightPos;
                void main() {
                   highp vec3 L = normalize(lightPos - vert);
                   highp float NL = max(dot(normalize(vertNormal), L), 0.0);
                   highp vec3 color = vec3(1.0, 0.0, 0.0);
                   highp vec3 col = clamp(color * 0.2 + color * 0.8 * NL, 0.0, 1.0);
                   fragColor = vec4(col, 1.0);
                }"""

    def vertex_shader_source(self):
        return """attribute vec4 vertex;
                attribute vec3 normal;
                varying vec3 vert;
                varying vec3 vertNormal;
                uniform mat4 projMatrix;
                uniform mat4 mvMatrix;
                uniform mat3 normalMatrix;
                void main() {
                   vert = vertex.xyz;
                   vertNormal = normalMatrix * normal;
                   gl_Position = projMatrix * mvMatrix * vertex;
                }"""

    def fragment_shader_source(self):
        return """varying highp vec3 vert;
                varying highp vec3 vertNormal;
                uniform highp vec3 lightPos;
                void main() {
                   highp vec3 L = normalize(lightPos - vert);
                   highp float NL = max(dot(normalize(vertNormal), L), 0.0);
                   highp vec3 color = vec3(0.39, 1.0, 0.0);
                   highp vec3 col = clamp(color * 0.2 + color * 0.8 * NL, 0.0, 1.0);
                   gl_FragColor = vec4(col, 1.0);
                }"""

    def initializeGL(self):
        self.context().aboutToBeDestroyed.connect(self.cleanup)
        self.initializeOpenGLFunctions()
        self.glClearColor(1.0, 1.0, 1.0, 0 if self._transparent else 1)

        self.program = QOpenGLShaderProgram()

        if self._core:
            self._vertex_shader = self.vertex_shader_source_core()
            self._fragment_shader = self.fragment_shader_source_core()
        else:
            self._vertex_shader = self.vertex_shader_source()
            self._fragment_shader = self.fragment_shader_source()

        self.program.addShaderFromSourceCode(QOpenGLShader.Vertex, self._vertex_shader)
        self.program.addShaderFromSourceCode(QOpenGLShader.Fragment, self._fragment_shader)
        self.program.bindAttributeLocation("vertex", 0)
        self.program.bindAttributeLocation("normal", 1)
        self.program.link()

        self.program.bind()
        self._proj_matrix_loc = self.program.uniformLocation("projMatrix")
        self._mv_matrix_loc = self.program.uniformLocation("mvMatrix")
        self._normal_matrix_loc = self.program.uniformLocation("normalMatrix")
        self._light_pos_loc = self.program.uniformLocation("lightPos")

        self.vao.create()
        with QOpenGLVertexArrayObject.Binder(self.vao):
            self._solver_vbo.create()
            self._solver_vbo.bind()
            float_size = ctypes.sizeof(ctypes.c_float)
            self._solver_vbo.allocate(self.fluid_solver.const_data(), self.fluid_solver.count() * float_size)

            self.setup_vertex_attribs()

            self.camera.setToIdentity()
            self.camera.translate(0, 0, -1)

            self.program.setUniformValue(self._light_pos_loc, QVector3D(0, 0, 70))
            self.program.release()

    def setup_vertex_attribs(self):
        self._solver_vbo.bind()
        f = QOpenGLContext.currentContext().functions()
        f.glEnableVertexAttribArray(0)
        f.glEnableVertexAttribArray(1)
        float_size = ctypes.sizeof(ctypes.c_float)

        null = VoidPtr(0)
        pointer = VoidPtr(3 * float_size)
        f.glVertexAttribPointer(0, 3, int(GL.GL_FLOAT), int(GL.GL_FALSE), 6 * float_size, null)
        f.glVertexAttribPointer(1, 3, int(GL.GL_FLOAT), int(GL.GL_FALSE), 6 * float_size, pointer)
        self._solver_vbo.release()

    def update_vbo(self):
        self._solver_vbo.bind()
        float_size = ctypes.sizeof(ctypes.c_float)
        self._solver_vbo.write(0, self.fluid_solver.const_data(), self.fluid_solver.count() * float_size)
        self._solver_vbo.release()
        
    def paintGL(self):

        self.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        self.glEnable(GL.GL_DEPTH_TEST)
        self.glEnable(GL.GL_CULL_FACE)

        self.world.setToIdentity()
        self.world.rotate(180 - (self._x_rot / 16), 1, 0, 0)
        self.world.rotate(self._y_rot / 16, 0, 1, 0)
        self.world.rotate(self._z_rot / 16, 0, 0, 1)
        self.world.scale(0.05 * self.zoom, 0.05 * self.zoom, 0.05 * self.zoom)
        self.world.translate(self._x_pan, self._y_pan, 0.0)
        
        self.update_vbo()

        if self.system_obj.start_playforward:
            self.current_frame_index = (self.current_frame_index + 1) % self.system_obj.num_frames
        elif self.system_obj.start_playback:
            self.current_frame_index = (self.current_frame_index - 1 + self.system_obj.num_frames) % self.system_obj.num_frames

        with QOpenGLVertexArrayObject.Binder(self.vao):
            self.program.bind()

            self.program.setUniformValue(self._proj_matrix_loc, self.proj)
            self.program.setUniformValue(self._mv_matrix_loc, self.camera * self.world)
            normal_matrix = self.world.normalMatrix()
            self.program.setUniformValue(self._normal_matrix_loc, normal_matrix)
            
            self.glDrawArrays(GL.GL_POINTS, 0, self.fluid_solver.count())

            self.program.release()

    def animate(self):
        if self.system_obj.stop == False:
            if self.system_obj.start_play == True and \
                self.system_obj.finished_caching == True:

                self.fluid_solver.update_m_data(self.current_frame_index)

                self.update()

    def resizeGL(self, width, height):
        self.proj.setToIdentity()
        self.proj.perspective(45, width / height, 0.01, 10000)

    def mousePressEvent(self, event):
        self._last_pos = event.position()

    def mouseMoveEvent(self, event):
        pos = event.position()
        dx = pos.x() - self._last_pos.x()
        dy = pos.y() - self._last_pos.y()

        if event.buttons() & Qt.MouseButton.LeftButton:
            self.set_xrotation(self._x_rot + 8 * dy)
            self.set_yrotation(self._y_rot + 8 * dx)
        elif event.buttons() & Qt.MouseButton.RightButton:
            self.set_xrotation(self._x_rot + 8 * dy)
            self.set_zrotation(self._z_rot + 8 * dx)
        elif event.buttons() & Qt.MouseButton.MiddleButton:
            self.set_xpan(self._x_pan + dx)
            self.set_ypan(self._y_pan + dy)

        self._last_pos = pos

    def wheelEvent(self, event: QWheelEvent) -> None:
        delta = event.angleDelta().y() / 120  # Normalized delta
        self.zoom += delta * 0.1  # Adjust the zoom speed

        # Ensure zoom stays within a reasonable range
        self.zoom = max(0.05, min(5.0, self.zoom))

        self.update()  # Trigger a repaint