# 
TARGET = libGEL.so

BIN_DIR=build
LIBSRCS = $(wildcard */*.cpp) $(wildcard */*.c)
LIBOBJS = $(LIBSRCS:.cpp=.o) 

	
CXX = g++
CXXFLAGS = -Wall -O3 -std=c++11
CFLAGS=-I ./ -I ../

$(TARGET): $(LIBOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(LIBOBJS)  $(CFLAGS) -L./Lib/libGLEW.1.5.dylib
	
%.o : %.cpp
	$(CXX) $(CXXFLAGS)  -c $< -o $@ $(CFLAGS)

.PHONY: clean test
clean: 
	rm $(LIBOBJS)
	
depend: $(LIBSRCS)
	makedepend -Y $(LIBSRCS) 2> /dev/null
# DO NOT DELETE

CGLA/ArithSqMat3x3Float.o: CGLA/ArithSqMat3x3Float.h CGLA/ExceptionStandard.h
CGLA/ArithSqMat3x3Float.o: CGLA/CGLA-util.h CGLA/Vec3f.h
CGLA/ArithSqMat3x3Float.o: CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h
CGLA/ArithSqMat3x3Float.o: CGLA/ArithVec.h CGLA/Vec3i.h CGLA/ArithVec3Int.h
CGLA/ArithSqMat3x3Float.o: CGLA/ArithVecInt.h CGLA/Vec3usi.h
CGLA/ArithSqMat3x3Float.o: CGLA/ArithSqMatFloat.h CGLA/ArithMatFloat.h
CGLA/ArithSqMat3x3Float.o: CGLA/Mat3x3f.h CGLA/Mat3x3d.h CGLA/Vec3d.h
CGLA/ArithSqMat4x4Float.o: CGLA/ArithSqMat4x4Float.h CGLA/ExceptionStandard.h
CGLA/ArithSqMat4x4Float.o: CGLA/CGLA-util.h CGLA/Vec3f.h
CGLA/ArithSqMat4x4Float.o: CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h
CGLA/ArithSqMat4x4Float.o: CGLA/ArithVec.h CGLA/Vec3i.h CGLA/ArithVec3Int.h
CGLA/ArithSqMat4x4Float.o: CGLA/ArithVecInt.h CGLA/Vec3usi.h CGLA/Vec4f.h
CGLA/ArithSqMat4x4Float.o: CGLA/ArithVec4Float.h CGLA/ArithSqMatFloat.h
CGLA/ArithSqMat4x4Float.o: CGLA/ArithMatFloat.h CGLA/Mat4x4f.h CGLA/Mat4x4d.h
CGLA/ArithSqMat4x4Float.o: CGLA/Vec3d.h CGLA/Vec4d.h
CGLA/ArithVec2Float.o: CGLA/Vec2f.h CGLA/ArithVec2Float.h
CGLA/ArithVec2Float.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h CGLA/CGLA-util.h
CGLA/ArithVec2Float.o: CGLA/Vec2i.h CGLA/Vec2d.h CGLA/Mat2x2f.h
CGLA/ArithVec2Float.o: CGLA/ArithSqMat2x2Float.h CGLA/ExceptionStandard.h
CGLA/ArithVec2Float.o: CGLA/ArithSqMatFloat.h CGLA/ArithMatFloat.h
CGLA/ArithVec2Float.o: CGLA/Mat2x2d.h
CGLA/ArithVec3Float.o: CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h
CGLA/ArithVec3Float.o: CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec3f.h
CGLA/ArithVec3Float.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
CGLA/ArithVec3Float.o: CGLA/Vec3usi.h CGLA/Vec3d.h
CGLA/Mat3x3d.o: CGLA/Mat3x3d.h CGLA/CGLA-util.h CGLA/Vec3d.h
CGLA/Mat3x3d.o: CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h CGLA/ArithVec.h
CGLA/Mat3x3d.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
CGLA/Mat3x3d.o: CGLA/Vec3usi.h CGLA/Vec3f.h CGLA/ArithSqMat3x3Float.h
CGLA/Mat3x3d.o: CGLA/ExceptionStandard.h CGLA/ArithSqMatFloat.h
CGLA/Mat3x3d.o: CGLA/ArithMatFloat.h
CGLA/Mat3x3f.o: CGLA/Mat3x3f.h CGLA/CGLA-util.h CGLA/Vec3f.h
CGLA/Mat3x3f.o: CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h CGLA/ArithVec.h
CGLA/Mat3x3f.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
CGLA/Mat3x3f.o: CGLA/Vec3usi.h CGLA/ArithSqMat3x3Float.h
CGLA/Mat3x3f.o: CGLA/ExceptionStandard.h CGLA/ArithSqMatFloat.h
CGLA/Mat3x3f.o: CGLA/ArithMatFloat.h
CGLA/Mat4x4d.o: CGLA/Mat4x4d.h CGLA/ExceptionStandard.h CGLA/CGLA-util.h
CGLA/Mat4x4d.o: CGLA/Vec3d.h CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h
CGLA/Mat4x4d.o: CGLA/ArithVec.h CGLA/Vec3i.h CGLA/ArithVec3Int.h
CGLA/Mat4x4d.o: CGLA/ArithVecInt.h CGLA/Vec3usi.h CGLA/Vec3f.h CGLA/Vec4d.h
CGLA/Mat4x4d.o: CGLA/ArithVec4Float.h CGLA/ArithSqMat4x4Float.h CGLA/Vec4f.h
CGLA/Mat4x4d.o: CGLA/ArithSqMatFloat.h CGLA/ArithMatFloat.h
CGLA/Mat4x4f.o: CGLA/Mat4x4f.h CGLA/ExceptionStandard.h CGLA/CGLA-util.h
CGLA/Mat4x4f.o: CGLA/Vec3f.h CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h
CGLA/Mat4x4f.o: CGLA/ArithVec.h CGLA/Vec3i.h CGLA/ArithVec3Int.h
CGLA/Mat4x4f.o: CGLA/ArithVecInt.h CGLA/Vec3usi.h CGLA/Vec4f.h
CGLA/Mat4x4f.o: CGLA/ArithVec4Float.h CGLA/ArithSqMat4x4Float.h
CGLA/Mat4x4f.o: CGLA/ArithSqMatFloat.h CGLA/ArithMatFloat.h
CGLA/TableTrigonometry.o: CGLA/TableTrigonometry.h CGLA/CGLA-util.h
CGLA/Vec2i.o: CGLA/Vec2i.h CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec2f.h
CGLA/Vec2i.o: CGLA/ArithVec2Float.h CGLA/ArithVecFloat.h
CGLA/Vec3f.o: CGLA/Vec3f.h CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h
CGLA/Vec3f.o: CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec3i.h
CGLA/Vec3f.o: CGLA/ArithVec3Int.h CGLA/ArithVecInt.h CGLA/Vec3usi.h
CGLA/Vec3f.o: CGLA/Vec3d.h CGLA/Vec4f.h CGLA/ArithVec4Float.h
CGLA/Vec3i.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
CGLA/Vec3i.o: CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec3f.h
CGLA/Vec3i.o: CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h CGLA/Vec3usi.h
CGLA/Vec3i.o: CGLA/Vec3d.h CGLA/Vec3uc.h
CGLA/eigensolution.o: CGLA/eigensolution.h CGLA/Mat2x2f.h CGLA/Vec2f.h
CGLA/eigensolution.o: CGLA/ArithVec2Float.h CGLA/ArithVecFloat.h
CGLA/eigensolution.o: CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec2i.h
CGLA/eigensolution.o: CGLA/ArithSqMat2x2Float.h CGLA/ExceptionStandard.h
CGLA/eigensolution.o: CGLA/ArithSqMatFloat.h CGLA/ArithMatFloat.h
CGLA/eigensolution.o: CGLA/Mat3x3f.h CGLA/Vec3f.h CGLA/ArithVec3Float.h
CGLA/eigensolution.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
CGLA/eigensolution.o: CGLA/Vec3usi.h CGLA/ArithSqMat3x3Float.h CGLA/Mat4x4f.h
CGLA/eigensolution.o: CGLA/Vec4f.h CGLA/ArithVec4Float.h
CGLA/eigensolution.o: CGLA/ArithSqMat4x4Float.h CGLA/Mat2x2d.h CGLA/Vec2d.h
CGLA/eigensolution.o: CGLA/Mat3x3d.h CGLA/Vec3d.h CGLA/Mat4x4d.h CGLA/Vec4d.h
CGLA/gel_rand.o: CGLA/CGLA-util.h
CGLA/statistics.o: CGLA/statistics.h CGLA/Mat2x2f.h CGLA/Vec2f.h
CGLA/statistics.o: CGLA/ArithVec2Float.h CGLA/ArithVecFloat.h CGLA/ArithVec.h
CGLA/statistics.o: CGLA/CGLA-util.h CGLA/Vec2i.h CGLA/ArithSqMat2x2Float.h
CGLA/statistics.o: CGLA/ExceptionStandard.h CGLA/ArithSqMatFloat.h
CGLA/statistics.o: CGLA/ArithMatFloat.h CGLA/Mat3x3f.h CGLA/Vec3f.h
CGLA/statistics.o: CGLA/ArithVec3Float.h CGLA/Vec3i.h CGLA/ArithVec3Int.h
CGLA/statistics.o: CGLA/ArithVecInt.h CGLA/Vec3usi.h
CGLA/statistics.o: CGLA/ArithSqMat3x3Float.h CGLA/Mat4x4f.h CGLA/Vec4f.h
CGLA/statistics.o: CGLA/ArithVec4Float.h CGLA/ArithSqMat4x4Float.h
CGLA/statistics.o: CGLA/Mat2x2d.h CGLA/Vec2d.h CGLA/Mat3x3d.h CGLA/Vec3d.h
CGLA/statistics.o: CGLA/Mat4x4d.h CGLA/Vec4d.h
GLGraphics/Console.o: GLGraphics/Console.h GL/glew.h GLGraphics/stb_image.h
GLGraphics/Console.o: GLGraphics/stb_image_write.h
GLGraphics/GLViewController.o: GL/glew.h GLGraphics/GLViewController.h
GLGraphics/GLViewController.o: GLGraphics/QuatTrackBall.h CGLA/Vec2i.h
GLGraphics/GLViewController.o: CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec2f.h
GLGraphics/GLViewController.o: CGLA/ArithVec2Float.h CGLA/ArithVecFloat.h
GLGraphics/GLViewController.o: CGLA/Vec2i.h CGLA/Vec3f.h
GLGraphics/GLViewController.o: CGLA/ArithVec3Float.h CGLA/Vec3i.h
GLGraphics/GLViewController.o: CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
GLGraphics/GLViewController.o: CGLA/Vec3usi.h CGLA/Quatf.h CGLA/ArithQuat.h
GLGraphics/GLViewController.o: CGLA/ArithSqMatFloat.h CGLA/ArithMatFloat.h
GLGraphics/GLViewController.o: CGLA/Vec3f.h CGLA/Vec4f.h
GLGraphics/GLViewController.o: CGLA/ArithVec4Float.h CGLA/Mat3x3f.h
GLGraphics/GLViewController.o: CGLA/ArithSqMat3x3Float.h
GLGraphics/GLViewController.o: CGLA/ExceptionStandard.h CGLA/Mat4x4f.h
GLGraphics/GLViewController.o: CGLA/ArithSqMat4x4Float.h CGLA/Mat3x3f.h
GLGraphics/IDBufferWireFrameRenderer.o: GLGraphics/IDBufferWireFrameRenderer.h
GLGraphics/IDBufferWireFrameRenderer.o: GL/glew.h CGLA/Vec3uc.h CGLA/Vec3i.h
GLGraphics/IDBufferWireFrameRenderer.o: CGLA/ArithVec3Int.h
GLGraphics/IDBufferWireFrameRenderer.o: CGLA/ArithVecInt.h CGLA/ArithVec.h
GLGraphics/IDBufferWireFrameRenderer.o: CGLA/CGLA-util.h CGLA/Vec3f.h
GLGraphics/IDBufferWireFrameRenderer.o: CGLA/ArithVec3Float.h
GLGraphics/IDBufferWireFrameRenderer.o: CGLA/ArithVecFloat.h CGLA/Vec3usi.h
GLGraphics/IDBufferWireFrameRenderer.o: HMesh/Manifold.h CGLA/Vec3d.h
GLGraphics/IDBufferWireFrameRenderer.o: CGLA/Vec3f.h
GLGraphics/IDBufferWireFrameRenderer.o: HMesh/ConnectivityKernel.h
GLGraphics/IDBufferWireFrameRenderer.o: HMesh/ItemVector.h HMesh/ItemID.h
GLGraphics/IDBufferWireFrameRenderer.o: HMesh/Iterators.h HMesh/Walker.h
GLGraphics/IDBufferWireFrameRenderer.o: HMesh/Manifold.h
GLGraphics/IDBufferWireFrameRenderer.o: HMesh/AttributeVector.h CGLA/Vec4f.h
GLGraphics/IDBufferWireFrameRenderer.o: CGLA/ArithVec4Float.h CGLA/Vec2f.h
GLGraphics/IDBufferWireFrameRenderer.o: CGLA/ArithVec2Float.h CGLA/Vec2i.h
GLGraphics/IDBufferWireFrameRenderer.o: HMesh/AttributeVector.h
GLGraphics/IDBufferWireFrameRenderer.o: GLGraphics/draw.h Geometry/TriMesh.h
GLGraphics/IDBufferWireFrameRenderer.o: CGLA/Mat4x4f.h
GLGraphics/IDBufferWireFrameRenderer.o: CGLA/ExceptionStandard.h CGLA/Vec4f.h
GLGraphics/IDBufferWireFrameRenderer.o: CGLA/ArithSqMat4x4Float.h
GLGraphics/IDBufferWireFrameRenderer.o: CGLA/ArithSqMatFloat.h
GLGraphics/IDBufferWireFrameRenderer.o: CGLA/ArithMatFloat.h
GLGraphics/IDBufferWireFrameRenderer.o: Geometry/IndexedFaceSet.h
GLGraphics/IDBufferWireFrameRenderer.o: CGLA/Vec3i.h Geometry/Material.h
GLGraphics/IDBufferWireFrameRenderer.o: Geometry/AABox.h Geometry/Triangle.h
GLGraphics/IDBufferWireFrameRenderer.o: CGLA/Mat3x3f.h Geometry/OBox.h
GLGraphics/IDBufferWireFrameRenderer.o: Geometry/AABox.h
GLGraphics/IDBufferWireFrameRenderer.o: Geometry/BoundingINode.h
GLGraphics/IDBufferWireFrameRenderer.o: Geometry/Ray.h Geometry/TriMesh.h
GLGraphics/IDBufferWireFrameRenderer.o: Geometry/BoundingNode.h
GLGraphics/IDBufferWireFrameRenderer.o: Geometry/OBox.h
GLGraphics/IDBufferWireFrameRenderer.o: Geometry/BoundingLNode.h
GLGraphics/IDBufferWireFrameRenderer.o: Geometry/BoundingTree.h
GLGraphics/IDBufferWireFrameRenderer.o: Geometry/BoundingLNode.h
GLGraphics/IDBufferWireFrameRenderer.o: Geometry/BoundingINode.h
GLGraphics/IDBufferWireFrameRenderer.o: Geometry/Graph.h Util/Range.h
GLGraphics/IDBufferWireFrameRenderer.o: Util/AttribVec.h
GLGraphics/IDBufferWireFrameRenderer.o: GLGraphics/glsl_shader.h
GLGraphics/ManifoldRenderer.o: GLGraphics/ManifoldRenderer.h CGLA/Vec2d.h
GLGraphics/ManifoldRenderer.o: CGLA/ArithVec2Float.h CGLA/ArithVecFloat.h
GLGraphics/ManifoldRenderer.o: CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec2i.h
GLGraphics/ManifoldRenderer.o: CGLA/Vec2f.h GL/glew.h GLGraphics/draw.h
GLGraphics/ManifoldRenderer.o: Geometry/TriMesh.h CGLA/Mat4x4f.h
GLGraphics/ManifoldRenderer.o: CGLA/ExceptionStandard.h CGLA/Vec3f.h
GLGraphics/ManifoldRenderer.o: CGLA/ArithVec3Float.h CGLA/Vec3i.h
GLGraphics/ManifoldRenderer.o: CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
GLGraphics/ManifoldRenderer.o: CGLA/Vec3usi.h CGLA/Vec4f.h
GLGraphics/ManifoldRenderer.o: CGLA/ArithVec4Float.h
GLGraphics/ManifoldRenderer.o: CGLA/ArithSqMat4x4Float.h
GLGraphics/ManifoldRenderer.o: CGLA/ArithSqMatFloat.h CGLA/ArithMatFloat.h
GLGraphics/ManifoldRenderer.o: Geometry/IndexedFaceSet.h CGLA/Vec3f.h
GLGraphics/ManifoldRenderer.o: CGLA/Vec3i.h Geometry/Material.h
GLGraphics/ManifoldRenderer.o: Geometry/AABox.h Geometry/Triangle.h
GLGraphics/ManifoldRenderer.o: CGLA/Vec2f.h CGLA/Mat3x3f.h Geometry/OBox.h
GLGraphics/ManifoldRenderer.o: Geometry/AABox.h Geometry/BoundingINode.h
GLGraphics/ManifoldRenderer.o: Geometry/Ray.h Geometry/TriMesh.h
GLGraphics/ManifoldRenderer.o: Geometry/BoundingNode.h Geometry/OBox.h
GLGraphics/ManifoldRenderer.o: Geometry/BoundingLNode.h
GLGraphics/ManifoldRenderer.o: Geometry/BoundingTree.h
GLGraphics/ManifoldRenderer.o: Geometry/BoundingLNode.h
GLGraphics/ManifoldRenderer.o: Geometry/BoundingINode.h Geometry/Graph.h
GLGraphics/ManifoldRenderer.o: CGLA/Vec3d.h Util/Range.h Util/AttribVec.h
GLGraphics/ManifoldRenderer.o: GLGraphics/Console.h
GLGraphics/ManifoldRenderer.o: GLGraphics/IDBufferWireFrameRenderer.h
GLGraphics/ManifoldRenderer.o: CGLA/Vec3uc.h HMesh/Manifold.h
GLGraphics/ManifoldRenderer.o: HMesh/ConnectivityKernel.h HMesh/ItemVector.h
GLGraphics/ManifoldRenderer.o: HMesh/ItemID.h HMesh/Iterators.h
GLGraphics/ManifoldRenderer.o: HMesh/Walker.h HMesh/Manifold.h
GLGraphics/ManifoldRenderer.o: HMesh/AttributeVector.h CGLA/Vec4d.h
GLGraphics/ManifoldRenderer.o: CGLA/Vec3d.h CGLA/Mat3x3d.h
GLGraphics/ManifoldRenderer.o: CGLA/ArithSqMat3x3Float.h
GLGraphics/ManifoldRenderer.o: GLGraphics/glsl_shader.h
GLGraphics/ManifoldRenderer.o: HMesh/AttributeVector.h HMesh/curvature.h
GLGraphics/MeshEditor.o: GLGraphics/MeshEditor.h GLGraphics/Console.h
GLGraphics/MeshEditor.o: GL/glew.h GLGraphics/VisObj.h HMesh/Manifold.h
GLGraphics/MeshEditor.o: CGLA/Vec3d.h CGLA/ArithVec3Float.h
GLGraphics/MeshEditor.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h
GLGraphics/MeshEditor.o: CGLA/CGLA-util.h CGLA/Vec3i.h CGLA/ArithVec3Int.h
GLGraphics/MeshEditor.o: CGLA/ArithVecInt.h CGLA/Vec3usi.h CGLA/Vec3f.h
GLGraphics/MeshEditor.o: HMesh/ConnectivityKernel.h HMesh/ItemVector.h
GLGraphics/MeshEditor.o: HMesh/ItemID.h HMesh/Iterators.h HMesh/Walker.h
GLGraphics/MeshEditor.o: HMesh/Manifold.h HMesh/AttributeVector.h
GLGraphics/MeshEditor.o: Geometry/Graph.h Util/Range.h Util/AttribVec.h
GLGraphics/MeshEditor.o: GLGraphics/draw.h Geometry/TriMesh.h CGLA/Mat4x4f.h
GLGraphics/MeshEditor.o: CGLA/ExceptionStandard.h CGLA/Vec4f.h
GLGraphics/MeshEditor.o: CGLA/ArithVec4Float.h CGLA/ArithSqMat4x4Float.h
GLGraphics/MeshEditor.o: CGLA/ArithSqMatFloat.h CGLA/ArithMatFloat.h
GLGraphics/MeshEditor.o: Geometry/IndexedFaceSet.h CGLA/Vec3f.h CGLA/Vec3i.h
GLGraphics/MeshEditor.o: Geometry/Material.h Geometry/AABox.h
GLGraphics/MeshEditor.o: Geometry/Triangle.h CGLA/Vec2f.h
GLGraphics/MeshEditor.o: CGLA/ArithVec2Float.h CGLA/Vec2i.h CGLA/Mat3x3f.h
GLGraphics/MeshEditor.o: Geometry/OBox.h Geometry/AABox.h
GLGraphics/MeshEditor.o: Geometry/BoundingINode.h Geometry/Ray.h
GLGraphics/MeshEditor.o: Geometry/TriMesh.h Geometry/BoundingNode.h
GLGraphics/MeshEditor.o: Geometry/OBox.h Geometry/BoundingLNode.h
GLGraphics/MeshEditor.o: Geometry/BoundingTree.h Geometry/BoundingLNode.h
GLGraphics/MeshEditor.o: Geometry/BoundingINode.h
GLGraphics/MeshEditor.o: GLGraphics/GLViewController.h
GLGraphics/MeshEditor.o: GLGraphics/QuatTrackBall.h CGLA/Vec2i.h CGLA/Quatf.h
GLGraphics/MeshEditor.o: CGLA/ArithQuat.h CGLA/Mat3x3f.h
GLGraphics/MeshEditor.o: CGLA/ArithSqMat3x3Float.h CGLA/Mat4x4f.h
GLGraphics/MeshEditor.o: GLGraphics/ManifoldRenderer.h CGLA/Vec2d.h
GLGraphics/MeshEditor.o: CGLA/Vec2f.h GLGraphics/IDBufferWireFrameRenderer.h
GLGraphics/MeshEditor.o: CGLA/Vec3uc.h CGLA/Vec4d.h CGLA/Vec3d.h
GLGraphics/MeshEditor.o: GLGraphics/VisObj.h
GLGraphics/QuatTrackBall.o: GL/glew.h CGLA/CGLA-util.h
GLGraphics/QuatTrackBall.o: GLGraphics/QuatTrackBall.h CGLA/Vec2i.h
GLGraphics/QuatTrackBall.o: CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec2f.h
GLGraphics/QuatTrackBall.o: CGLA/ArithVec2Float.h CGLA/ArithVecFloat.h
GLGraphics/QuatTrackBall.o: CGLA/Vec2i.h CGLA/Vec3f.h CGLA/ArithVec3Float.h
GLGraphics/QuatTrackBall.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h
GLGraphics/QuatTrackBall.o: CGLA/ArithVecInt.h CGLA/Vec3usi.h CGLA/Quatf.h
GLGraphics/QuatTrackBall.o: CGLA/ArithQuat.h CGLA/ArithSqMatFloat.h
GLGraphics/QuatTrackBall.o: CGLA/ArithMatFloat.h CGLA/Vec3f.h CGLA/Vec4f.h
GLGraphics/QuatTrackBall.o: CGLA/ArithVec4Float.h CGLA/Mat3x3f.h
GLGraphics/QuatTrackBall.o: CGLA/ArithSqMat3x3Float.h
GLGraphics/QuatTrackBall.o: CGLA/ExceptionStandard.h CGLA/Mat4x4f.h
GLGraphics/QuatTrackBall.o: CGLA/ArithSqMat4x4Float.h
GLGraphics/QuatTrackBall.o: GLGraphics/GLViewController.h
GLGraphics/ShadowBuffer.o: GLGraphics/ShadowBuffer.h
GLGraphics/SimpleTrackBall.o: GL/glew.h GLGraphics/GLViewController.h
GLGraphics/SimpleTrackBall.o: GLGraphics/QuatTrackBall.h CGLA/Vec2i.h
GLGraphics/SimpleTrackBall.o: CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec2f.h
GLGraphics/SimpleTrackBall.o: CGLA/ArithVec2Float.h CGLA/ArithVecFloat.h
GLGraphics/SimpleTrackBall.o: CGLA/Vec2i.h CGLA/Vec3f.h CGLA/ArithVec3Float.h
GLGraphics/SimpleTrackBall.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h
GLGraphics/SimpleTrackBall.o: CGLA/ArithVecInt.h CGLA/Vec3usi.h CGLA/Quatf.h
GLGraphics/SimpleTrackBall.o: CGLA/ArithQuat.h CGLA/ArithSqMatFloat.h
GLGraphics/SimpleTrackBall.o: CGLA/ArithMatFloat.h CGLA/Vec3f.h CGLA/Vec4f.h
GLGraphics/SimpleTrackBall.o: CGLA/ArithVec4Float.h CGLA/Mat3x3f.h
GLGraphics/SimpleTrackBall.o: CGLA/ArithSqMat3x3Float.h
GLGraphics/SimpleTrackBall.o: CGLA/ExceptionStandard.h CGLA/Mat4x4f.h
GLGraphics/SimpleTrackBall.o: CGLA/ArithSqMat4x4Float.h
GLGraphics/SimpleTrackBall.o: GLGraphics/SimpleTrackBall.h CGLA/CGLA-util.h
GLGraphics/SinglePassWireframeRenderer.o: GLGraphics/glsl_shader.h GL/glew.h
GLGraphics/SinglePassWireframeRenderer.o: GLGraphics/SinglePassWireframeRenderer.h
GLGraphics/SinglePassWireframeRenderer.o: CGLA/Vec3f.h CGLA/ArithVec3Float.h
GLGraphics/SinglePassWireframeRenderer.o: CGLA/ArithVecFloat.h
GLGraphics/SinglePassWireframeRenderer.o: CGLA/ArithVec.h CGLA/CGLA-util.h
GLGraphics/SinglePassWireframeRenderer.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h
GLGraphics/SinglePassWireframeRenderer.o: CGLA/ArithVecInt.h CGLA/Vec3usi.h
GLGraphics/VisObj.o: GLGraphics/VisObj.h
GLGraphics/draw.o: GL/glew.h CGLA/Mat4x4f.h CGLA/ExceptionStandard.h
GLGraphics/draw.o: CGLA/CGLA-util.h CGLA/Vec3f.h CGLA/ArithVec3Float.h
GLGraphics/draw.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h CGLA/Vec3i.h
GLGraphics/draw.o: CGLA/ArithVec3Int.h CGLA/ArithVecInt.h CGLA/Vec3usi.h
GLGraphics/draw.o: CGLA/Vec4f.h CGLA/ArithVec4Float.h
GLGraphics/draw.o: CGLA/ArithSqMat4x4Float.h CGLA/ArithSqMatFloat.h
GLGraphics/draw.o: CGLA/ArithMatFloat.h CGLA/Vec3d.h HMesh/Manifold.h
GLGraphics/draw.o: HMesh/ConnectivityKernel.h HMesh/ItemVector.h
GLGraphics/draw.o: HMesh/ItemID.h HMesh/Iterators.h HMesh/Walker.h
GLGraphics/draw.o: HMesh/Manifold.h HMesh/AttributeVector.h
GLGraphics/draw.o: Geometry/TriMesh.h Geometry/IndexedFaceSet.h CGLA/Vec3f.h
GLGraphics/draw.o: CGLA/Vec3i.h Geometry/Material.h GLGraphics/draw.h
GLGraphics/draw.o: Geometry/AABox.h Geometry/Triangle.h CGLA/Vec2f.h
GLGraphics/draw.o: CGLA/ArithVec2Float.h CGLA/Vec2i.h CGLA/Mat3x3f.h
GLGraphics/draw.o: Geometry/OBox.h Geometry/AABox.h Geometry/BoundingINode.h
GLGraphics/draw.o: Geometry/Ray.h Geometry/TriMesh.h Geometry/BoundingNode.h
GLGraphics/draw.o: Geometry/OBox.h Geometry/BoundingLNode.h
GLGraphics/draw.o: Geometry/BoundingTree.h Geometry/BoundingLNode.h
GLGraphics/draw.o: Geometry/BoundingINode.h Geometry/Graph.h Util/Range.h
GLGraphics/draw.o: Util/AttribVec.h GLGraphics/SinglePassWireframeRenderer.h
GLGraphics/draw.o: GLGraphics/IDBufferWireFrameRenderer.h CGLA/Vec3uc.h
GLGraphics/draw.o: GLGraphics/SOIL.h
GLGraphics/glsl_shader.o: GL/glew.h
Geometry/AABox.o: Geometry/AABox.h
Geometry/BBox.o: Geometry/Ray.h CGLA/Vec3i.h CGLA/Vec3f.h
Geometry/BBox.o: CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h CGLA/ArithVec.h
Geometry/BBox.o: CGLA/CGLA-util.h CGLA/Vec3i.h CGLA/ArithVec3Int.h
Geometry/BBox.o: CGLA/ArithVecInt.h CGLA/Vec3usi.h Geometry/TriMesh.h
Geometry/BBox.o: Geometry/Material.h Geometry/BBox.h
Geometry/BSPTree.o: CGLA/Vec3d.h CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h
Geometry/BSPTree.o: CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec3i.h
Geometry/BSPTree.o: CGLA/ArithVec3Int.h CGLA/ArithVecInt.h CGLA/Vec3usi.h
Geometry/BSPTree.o: CGLA/Vec3f.h Geometry/TriMesh.h CGLA/Mat4x4f.h
Geometry/BSPTree.o: CGLA/ExceptionStandard.h CGLA/Vec4f.h
Geometry/BSPTree.o: CGLA/ArithVec4Float.h CGLA/ArithSqMat4x4Float.h
Geometry/BSPTree.o: CGLA/ArithSqMatFloat.h CGLA/ArithMatFloat.h
Geometry/BSPTree.o: Geometry/IndexedFaceSet.h CGLA/Vec3f.h CGLA/Vec3i.h
Geometry/BSPTree.o: Geometry/Material.h Geometry/BSPTree.h Geometry/BBox.h
Geometry/BSPTree.o: Geometry/Ray.h Geometry/TriMesh.h
Geometry/BoundingINode.o: Geometry/AABox.h Geometry/BoundingINode.h
Geometry/BoundingLNode.o: Geometry/AABox.h Geometry/BoundingLNode.h
Geometry/BoundingNode.o: Geometry/BoundingNode.h CGLA/Vec3f.h
Geometry/BoundingNode.o: CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h
Geometry/BoundingNode.o: CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec3i.h
Geometry/BoundingNode.o: CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
Geometry/BoundingNode.o: CGLA/Vec3usi.h Geometry/Ray.h CGLA/Vec3i.h
Geometry/BoundingNode.o: Geometry/TriMesh.h Geometry/Material.h
Geometry/BoundingNode.o: Geometry/Triangle.h CGLA/Vec2f.h
Geometry/BoundingNode.o: CGLA/ArithVec2Float.h CGLA/Vec2i.h CGLA/Mat3x3f.h
Geometry/BoundingNode.o: Geometry/AABox.h Geometry/OBox.h
Geometry/BoundingNode.o: Geometry/BoundingINode.h Geometry/BoundingLNode.h
Geometry/BoundingTree.o: Geometry/AABox.h Geometry/OBox.h
Geometry/BoundingTree.o: Geometry/BoundingTree.h Geometry/BoundingNode.h
Geometry/BoundingTree.o: CGLA/Vec3f.h CGLA/ArithVec3Float.h
Geometry/BoundingTree.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h
Geometry/BoundingTree.o: CGLA/CGLA-util.h CGLA/Vec3i.h CGLA/ArithVec3Int.h
Geometry/BoundingTree.o: CGLA/ArithVecInt.h CGLA/Vec3usi.h Geometry/Ray.h
Geometry/BoundingTree.o: CGLA/Vec3i.h Geometry/TriMesh.h Geometry/Material.h
Geometry/BoundingTree.o: Geometry/Triangle.h CGLA/Vec2f.h
Geometry/BoundingTree.o: CGLA/ArithVec2Float.h CGLA/Vec2i.h CGLA/Mat3x3f.h
Geometry/BoundingTree.o: Geometry/BoundingLNode.h Geometry/BoundingINode.h
Geometry/GradientFilter.o: Geometry/GradientFilter.h CGLA/Vec3i.h
Geometry/GradientFilter.o: CGLA/Vec3f.h CGLA/ArithVec3Float.h
Geometry/GradientFilter.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h
Geometry/GradientFilter.o: CGLA/CGLA-util.h CGLA/Vec3i.h CGLA/ArithVec3Int.h
Geometry/GradientFilter.o: CGLA/ArithVecInt.h CGLA/Vec3usi.h
Geometry/GradientFilter.o: Geometry/Neighbours.h CGLA/Vec3d.h CGLA/Vec3f.h
Geometry/GradientFilter.o: Geometry/RGrid.h Geometry/AncestorGrid.h
Geometry/Graph.o: Geometry/Graph.h CGLA/Vec3d.h CGLA/ArithVec3Float.h
Geometry/Graph.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h CGLA/CGLA-util.h
Geometry/Graph.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
Geometry/Graph.o: CGLA/Vec3usi.h CGLA/Vec3f.h Util/Range.h Util/AttribVec.h
Geometry/Implicit.o: CGLA/CGLA-util.h Geometry/Neighbours.h CGLA/Vec3i.h
Geometry/Implicit.o: CGLA/Vec3f.h CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h
Geometry/Implicit.o: CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec3i.h
Geometry/Implicit.o: CGLA/ArithVec3Int.h CGLA/ArithVecInt.h CGLA/Vec3usi.h
Geometry/Implicit.o: CGLA/Vec3d.h CGLA/Vec3f.h Geometry/GridAlgorithm.h
Geometry/Implicit.o: Geometry/RGrid.h Geometry/AncestorGrid.h
Geometry/Implicit.o: Geometry/HGrid.h Geometry/Cell.h CGLA/BitMask.h
Geometry/Implicit.o: Geometry/XForm.h Geometry/Implicit.h
Geometry/Neighbours.o: Geometry/Neighbours.h CGLA/Vec3i.h CGLA/Vec3f.h
Geometry/Neighbours.o: CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h
Geometry/Neighbours.o: CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec3i.h
Geometry/Neighbours.o: CGLA/ArithVec3Int.h CGLA/ArithVecInt.h CGLA/Vec3usi.h
Geometry/Neighbours.o: CGLA/Vec3d.h CGLA/Vec3f.h
Geometry/OBox.o: CGLA/statistics.h CGLA/eigensolution.h CGLA/Mat4x4f.h
Geometry/OBox.o: CGLA/ExceptionStandard.h CGLA/CGLA-util.h CGLA/Vec3f.h
Geometry/OBox.o: CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h CGLA/ArithVec.h
Geometry/OBox.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
Geometry/OBox.o: CGLA/Vec3usi.h CGLA/Vec4f.h CGLA/ArithVec4Float.h
Geometry/OBox.o: CGLA/ArithSqMat4x4Float.h CGLA/ArithSqMatFloat.h
Geometry/OBox.o: CGLA/ArithMatFloat.h Geometry/AABox.h Geometry/OBox.h
Geometry/OBox.o: Geometry/Triangle.h CGLA/Vec2f.h CGLA/ArithVec2Float.h
Geometry/OBox.o: CGLA/Vec2i.h CGLA/Vec3f.h CGLA/Mat3x3f.h
Geometry/Polygonizer.o: CGLA/CGLA-util.h Geometry/Polygonizer.h
Geometry/QEM.o: Geometry/QEM.h CGLA/Vec3d.h CGLA/ArithVec3Float.h
Geometry/QEM.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h CGLA/CGLA-util.h
Geometry/QEM.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
Geometry/QEM.o: CGLA/Vec3usi.h CGLA/Vec3f.h CGLA/Mat3x3d.h CGLA/Vec3d.h
Geometry/QEM.o: CGLA/ArithSqMat3x3Float.h CGLA/ExceptionStandard.h
Geometry/QEM.o: CGLA/ArithSqMatFloat.h CGLA/ArithMatFloat.h
Geometry/ThreeDDDA.o: CGLA/Vec3d.h CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h
Geometry/ThreeDDDA.o: CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec3i.h
Geometry/ThreeDDDA.o: CGLA/ArithVec3Int.h CGLA/ArithVecInt.h CGLA/Vec3usi.h
Geometry/ThreeDDDA.o: CGLA/Vec3f.h Geometry/ThreeDDDA.h CGLA/Vec3f.h
Geometry/ThreeDDDA.o: CGLA/Vec3i.h
Geometry/TriMesh.o: CGLA/Vec3f.h CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h
Geometry/TriMesh.o: CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec3i.h
Geometry/TriMesh.o: CGLA/ArithVec3Int.h CGLA/ArithVecInt.h CGLA/Vec3usi.h
Geometry/TriMesh.o: Geometry/TriMesh.h
Geometry/Triangle.o: CGLA/Vec4f.h CGLA/Vec3f.h CGLA/ArithVec3Float.h
Geometry/Triangle.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h CGLA/CGLA-util.h
Geometry/Triangle.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
Geometry/Triangle.o: CGLA/Vec3usi.h CGLA/ArithVec4Float.h Geometry/Triangle.h
Geometry/Triangle.o: CGLA/Vec2f.h CGLA/ArithVec2Float.h CGLA/Vec2i.h
Geometry/Triangle.o: CGLA/Vec3f.h CGLA/Mat3x3f.h
Geometry/TrilinFilter.o: Geometry/GradientFilter.h CGLA/Vec3i.h CGLA/Vec3f.h
Geometry/TrilinFilter.o: CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h
Geometry/TrilinFilter.o: CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec3i.h
Geometry/TrilinFilter.o: CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
Geometry/TrilinFilter.o: CGLA/Vec3usi.h Geometry/TrilinFilter.h
Geometry/TrilinFilter.o: Geometry/RGrid.h Geometry/AncestorGrid.h
Geometry/TrilinFilter.o: Geometry/Neighbours.h CGLA/Vec3d.h CGLA/Vec3f.h
Geometry/build_bbtree.o: Geometry/build_bbtree.h Geometry/BoundingTree.h
Geometry/build_bbtree.o: Geometry/BoundingNode.h CGLA/Vec3f.h
Geometry/build_bbtree.o: CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h
Geometry/build_bbtree.o: CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec3i.h
Geometry/build_bbtree.o: CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
Geometry/build_bbtree.o: CGLA/Vec3usi.h Geometry/Ray.h CGLA/Vec3i.h
Geometry/build_bbtree.o: Geometry/TriMesh.h Geometry/Material.h
Geometry/build_bbtree.o: Geometry/Triangle.h CGLA/Vec2f.h
Geometry/build_bbtree.o: CGLA/ArithVec2Float.h CGLA/Vec2i.h CGLA/Mat3x3f.h
Geometry/build_bbtree.o: Geometry/AABox.h Geometry/OBox.h
Geometry/build_bbtree.o: Geometry/BoundingLNode.h Geometry/BoundingINode.h
Geometry/build_bbtree.o: HMesh/Manifold.h CGLA/Vec3d.h CGLA/Vec3f.h
Geometry/build_bbtree.o: HMesh/ConnectivityKernel.h HMesh/ItemVector.h
Geometry/build_bbtree.o: HMesh/ItemID.h HMesh/Iterators.h HMesh/Walker.h
Geometry/build_bbtree.o: HMesh/Manifold.h HMesh/AttributeVector.h
Geometry/geo_load.o: Geometry/load.h Geometry/TriMesh.h Geometry/obj_load.h
Geometry/geo_load.o: Geometry/ply_load.h Geometry/TriMesh.h CGLA/Mat4x4f.h
Geometry/geo_load.o: CGLA/ExceptionStandard.h CGLA/CGLA-util.h CGLA/Vec3f.h
Geometry/geo_load.o: CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h
Geometry/geo_load.o: CGLA/ArithVec.h CGLA/Vec3i.h CGLA/ArithVec3Int.h
Geometry/geo_load.o: CGLA/ArithVecInt.h CGLA/Vec3usi.h CGLA/Vec4f.h
Geometry/geo_load.o: CGLA/ArithVec4Float.h CGLA/ArithSqMat4x4Float.h
Geometry/geo_load.o: CGLA/ArithSqMatFloat.h CGLA/ArithMatFloat.h
Geometry/geo_load.o: Geometry/IndexedFaceSet.h CGLA/Vec3f.h CGLA/Vec3i.h
Geometry/geo_load.o: Geometry/Material.h
Geometry/geo_obj_load.o: Geometry/TriMesh.h CGLA/Vec3f.h
Geometry/geo_obj_load.o: CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h
Geometry/geo_obj_load.o: CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec3i.h
Geometry/geo_obj_load.o: CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
Geometry/geo_obj_load.o: CGLA/Vec3usi.h
Geometry/load_raw.o: Geometry/load_raw.h CGLA/Mat4x4f.h
Geometry/load_raw.o: CGLA/ExceptionStandard.h CGLA/CGLA-util.h CGLA/Vec3f.h
Geometry/load_raw.o: CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h
Geometry/load_raw.o: CGLA/ArithVec.h CGLA/Vec3i.h CGLA/ArithVec3Int.h
Geometry/load_raw.o: CGLA/ArithVecInt.h CGLA/Vec3usi.h CGLA/Vec4f.h
Geometry/load_raw.o: CGLA/ArithVec4Float.h CGLA/ArithSqMat4x4Float.h
Geometry/load_raw.o: CGLA/ArithSqMatFloat.h CGLA/ArithMatFloat.h
Geometry/load_raw.o: Geometry/RGrid.h CGLA/Vec3i.h Geometry/AncestorGrid.h
Geometry/rply_load.o: Geometry/rply.h Geometry/ply_load.h Geometry/TriMesh.h
Geometry/rply_load.o: CGLA/Mat4x4f.h CGLA/ExceptionStandard.h
Geometry/rply_load.o: CGLA/CGLA-util.h CGLA/Vec3f.h CGLA/ArithVec3Float.h
Geometry/rply_load.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h CGLA/Vec3i.h
Geometry/rply_load.o: CGLA/ArithVec3Int.h CGLA/ArithVecInt.h CGLA/Vec3usi.h
Geometry/rply_load.o: CGLA/Vec4f.h CGLA/ArithVec4Float.h
Geometry/rply_load.o: CGLA/ArithSqMat4x4Float.h CGLA/ArithSqMatFloat.h
Geometry/rply_load.o: CGLA/ArithMatFloat.h Geometry/IndexedFaceSet.h
Geometry/rply_load.o: CGLA/Vec3f.h CGLA/Vec3i.h Geometry/Material.h
Geometry/rply_load.o: CGLA/Vec4f.h
Geometry/tessellate.o: Util/Grid2D.h CGLA/Vec2i.h CGLA/ArithVec.h
Geometry/tessellate.o: CGLA/CGLA-util.h CGLA/Vec3f.h CGLA/ArithVec3Float.h
Geometry/tessellate.o: CGLA/ArithVecFloat.h CGLA/Vec3i.h CGLA/ArithVec3Int.h
Geometry/tessellate.o: CGLA/ArithVecInt.h CGLA/Vec3usi.h CGLA/Vec2f.h
Geometry/tessellate.o: CGLA/ArithVec2Float.h CGLA/Vec2i.h
Geometry/tessellate.o: Geometry/tessellate.h Geometry/IndexedFaceSet.h
Geometry/tessellate.o: CGLA/Vec3i.h
Geometry/verification.o: Geometry/verification.h Geometry/Triangle.h
Geometry/verification.o: CGLA/Vec2f.h CGLA/ArithVec2Float.h
Geometry/verification.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h
Geometry/verification.o: CGLA/CGLA-util.h CGLA/Vec2i.h CGLA/Vec3f.h
Geometry/verification.o: CGLA/ArithVec3Float.h CGLA/Vec3i.h
Geometry/verification.o: CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
Geometry/verification.o: CGLA/Vec3usi.h CGLA/Mat3x3f.h
HMesh/ConnectivityKernel.o: HMesh/ConnectivityKernel.h HMesh/ItemVector.h
HMesh/ConnectivityKernel.o: HMesh/ItemID.h HMesh/Iterators.h
HMesh/Manifold.o: HMesh/Manifold.h CGLA/Vec3d.h CGLA/ArithVec3Float.h
HMesh/Manifold.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h CGLA/CGLA-util.h
HMesh/Manifold.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
HMesh/Manifold.o: CGLA/Vec3usi.h CGLA/Vec3f.h HMesh/ConnectivityKernel.h
HMesh/Manifold.o: HMesh/ItemVector.h HMesh/ItemID.h HMesh/Iterators.h
HMesh/Manifold.o: HMesh/Walker.h HMesh/AttributeVector.h Geometry/TriMesh.h
HMesh/Manifold.o: CGLA/Mat4x4f.h CGLA/ExceptionStandard.h CGLA/Vec4f.h
HMesh/Manifold.o: CGLA/ArithVec4Float.h CGLA/ArithSqMat4x4Float.h
HMesh/Manifold.o: CGLA/ArithSqMatFloat.h CGLA/ArithMatFloat.h
HMesh/Manifold.o: Geometry/IndexedFaceSet.h CGLA/Vec3f.h CGLA/Vec3i.h
HMesh/Manifold.o: Geometry/Material.h
HMesh/cleanup.o: HMesh/cleanup.h CGLA/Vec3f.h CGLA/ArithVec3Float.h
HMesh/cleanup.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h CGLA/CGLA-util.h
HMesh/cleanup.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
HMesh/cleanup.o: CGLA/Vec3usi.h CGLA/Vec3d.h CGLA/Vec3f.h Geometry/QEM.h
HMesh/cleanup.o: CGLA/Mat3x3d.h CGLA/Vec3d.h CGLA/ArithSqMat3x3Float.h
HMesh/cleanup.o: CGLA/ExceptionStandard.h CGLA/ArithSqMatFloat.h
HMesh/cleanup.o: CGLA/ArithMatFloat.h Geometry/KDTree.h CGLA/CGLA-util.h
HMesh/cleanup.o: CGLA/ArithVec.h HMesh/Manifold.h HMesh/ConnectivityKernel.h
HMesh/cleanup.o: HMesh/ItemVector.h HMesh/ItemID.h HMesh/Iterators.h
HMesh/cleanup.o: HMesh/Walker.h HMesh/AttributeVector.h
HMesh/curvature.o: HMesh/curvature.h HMesh/Manifold.h CGLA/Vec3d.h
HMesh/curvature.o: CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h CGLA/ArithVec.h
HMesh/curvature.o: CGLA/CGLA-util.h CGLA/Vec3i.h CGLA/ArithVec3Int.h
HMesh/curvature.o: CGLA/ArithVecInt.h CGLA/Vec3usi.h CGLA/Vec3f.h
HMesh/curvature.o: HMesh/ConnectivityKernel.h HMesh/ItemVector.h
HMesh/curvature.o: HMesh/ItemID.h HMesh/Iterators.h HMesh/Walker.h
HMesh/curvature.o: HMesh/AttributeVector.h CGLA/CGLA.h CGLA/BitMask.h
HMesh/curvature.o: CGLA/ExceptionStandard.h CGLA/Mat2x2d.h CGLA/Vec2d.h
HMesh/curvature.o: CGLA/ArithVec2Float.h CGLA/Vec2i.h CGLA/Vec2f.h
HMesh/curvature.o: CGLA/ArithSqMat2x2Float.h CGLA/ArithSqMatFloat.h
HMesh/curvature.o: CGLA/ArithMatFloat.h CGLA/Mat2x2f.h CGLA/Mat2x3d.h
HMesh/curvature.o: CGLA/Vec3d.h CGLA/Mat2x3f.h CGLA/Mat3x3d.h
HMesh/curvature.o: CGLA/ArithSqMat3x3Float.h CGLA/Mat3x3f.h CGLA/Mat4x4d.h
HMesh/curvature.o: CGLA/Vec4d.h CGLA/ArithVec4Float.h
HMesh/curvature.o: CGLA/ArithSqMat4x4Float.h CGLA/Vec4f.h CGLA/Mat4x4f.h
HMesh/curvature.o: CGLA/Quatd.h CGLA/ArithQuat.h CGLA/Quaternion.h
HMesh/curvature.o: CGLA/Quatf.h CGLA/TableTrigonometry.h CGLA/UnitVector.h
HMesh/curvature.o: CGLA/Vec2ui.h CGLA/Vec3uc.h CGLA/Vec4i.h
HMesh/curvature.o: CGLA/ArithVec4Int.h CGLA/Vec4uc.h CGLA/eigensolution.h
HMesh/curvature.o: CGLA/statistics.h CGLA/ls_solve.h HMesh/x3d_save.h
HMesh/curvature.o: HMesh/x3d_load.h Geometry/obj_load.h Geometry/TriMesh.h
HMesh/curvature.o: HMesh/mesh_optimization.h Geometry/Implicit.h
HMesh/curvature.o: Geometry/RGrid.h CGLA/Vec3i.h Geometry/AncestorGrid.h
HMesh/curvature.o: Geometry/XForm.h
HMesh/dual.o: HMesh/dual.h CGLA/Vec3d.h CGLA/ArithVec3Float.h
HMesh/dual.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h CGLA/CGLA-util.h
HMesh/dual.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
HMesh/dual.o: CGLA/Vec3usi.h CGLA/Vec3f.h Geometry/load.h Geometry/TriMesh.h
HMesh/dual.o: HMesh/Manifold.h HMesh/ConnectivityKernel.h HMesh/ItemVector.h
HMesh/dual.o: HMesh/ItemID.h HMesh/Iterators.h HMesh/Walker.h
HMesh/dual.o: HMesh/AttributeVector.h
HMesh/graph_algorithm.o: HMesh/graph_algorithm.h HMesh/Manifold.h
HMesh/graph_algorithm.o: CGLA/Vec3d.h CGLA/ArithVec3Float.h
HMesh/graph_algorithm.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h
HMesh/graph_algorithm.o: CGLA/CGLA-util.h CGLA/Vec3i.h CGLA/ArithVec3Int.h
HMesh/graph_algorithm.o: CGLA/ArithVecInt.h CGLA/Vec3usi.h CGLA/Vec3f.h
HMesh/graph_algorithm.o: HMesh/ConnectivityKernel.h HMesh/ItemVector.h
HMesh/graph_algorithm.o: HMesh/ItemID.h HMesh/Iterators.h HMesh/Walker.h
HMesh/graph_algorithm.o: HMesh/AttributeVector.h
HMesh/load.o: Geometry/load.h Geometry/TriMesh.h HMesh/Manifold.h
HMesh/load.o: CGLA/Vec3d.h CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h
HMesh/load.o: CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec3i.h
HMesh/load.o: CGLA/ArithVec3Int.h CGLA/ArithVecInt.h CGLA/Vec3usi.h
HMesh/load.o: CGLA/Vec3f.h HMesh/ConnectivityKernel.h HMesh/ItemVector.h
HMesh/load.o: HMesh/ItemID.h HMesh/Iterators.h HMesh/Walker.h
HMesh/load.o: HMesh/AttributeVector.h Geometry/ply_load.h Geometry/TriMesh.h
HMesh/load.o: CGLA/Mat4x4f.h CGLA/ExceptionStandard.h CGLA/Vec4f.h
HMesh/load.o: CGLA/ArithVec4Float.h CGLA/ArithSqMat4x4Float.h
HMesh/load.o: CGLA/ArithSqMatFloat.h CGLA/ArithMatFloat.h
HMesh/load.o: Geometry/IndexedFaceSet.h CGLA/Vec3f.h CGLA/Vec3i.h
HMesh/load.o: Geometry/Material.h HMesh/x3d_load.h Geometry/obj_load.h
HMesh/load.o: HMesh/off_load.h HMesh/cleanup.h
HMesh/mesh_optimization.o: HMesh/mesh_optimization.h HMesh/Manifold.h
HMesh/mesh_optimization.o: CGLA/Vec3d.h CGLA/ArithVec3Float.h
HMesh/mesh_optimization.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h
HMesh/mesh_optimization.o: CGLA/CGLA-util.h CGLA/Vec3i.h CGLA/ArithVec3Int.h
HMesh/mesh_optimization.o: CGLA/ArithVecInt.h CGLA/Vec3usi.h CGLA/Vec3f.h
HMesh/mesh_optimization.o: HMesh/ConnectivityKernel.h HMesh/ItemVector.h
HMesh/mesh_optimization.o: HMesh/ItemID.h HMesh/Iterators.h HMesh/Walker.h
HMesh/mesh_optimization.o: HMesh/AttributeVector.h Geometry/Implicit.h
HMesh/mesh_optimization.o: Geometry/RGrid.h CGLA/Vec3i.h
HMesh/mesh_optimization.o: Geometry/AncestorGrid.h Geometry/XForm.h
HMesh/mesh_optimization.o: HMesh/triangulate.h HMesh/smooth.h
HMesh/obj_load.o: Geometry/load.h Geometry/TriMesh.h Geometry/obj_load.h
HMesh/obj_load.o: HMesh/Manifold.h CGLA/Vec3d.h CGLA/ArithVec3Float.h
HMesh/obj_load.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h CGLA/CGLA-util.h
HMesh/obj_load.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
HMesh/obj_load.o: CGLA/Vec3usi.h CGLA/Vec3f.h HMesh/ConnectivityKernel.h
HMesh/obj_load.o: HMesh/ItemVector.h HMesh/ItemID.h HMesh/Iterators.h
HMesh/obj_load.o: HMesh/Walker.h HMesh/AttributeVector.h HMesh/cleanup.h
HMesh/obj_save.o: HMesh/obj_save.h CGLA/Vec3f.h CGLA/ArithVec3Float.h
HMesh/obj_save.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h CGLA/CGLA-util.h
HMesh/obj_save.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
HMesh/obj_save.o: CGLA/Vec3usi.h HMesh/Manifold.h CGLA/Vec3d.h CGLA/Vec3f.h
HMesh/obj_save.o: HMesh/ConnectivityKernel.h HMesh/ItemVector.h
HMesh/obj_save.o: HMesh/ItemID.h HMesh/Iterators.h HMesh/Walker.h
HMesh/obj_save.o: HMesh/AttributeVector.h
HMesh/off_load.o: HMesh/off_load.h Geometry/load.h Geometry/TriMesh.h
HMesh/off_load.o: HMesh/Manifold.h CGLA/Vec3d.h CGLA/ArithVec3Float.h
HMesh/off_load.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h CGLA/CGLA-util.h
HMesh/off_load.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
HMesh/off_load.o: CGLA/Vec3usi.h CGLA/Vec3f.h HMesh/ConnectivityKernel.h
HMesh/off_load.o: HMesh/ItemVector.h HMesh/ItemID.h HMesh/Iterators.h
HMesh/off_load.o: HMesh/Walker.h HMesh/AttributeVector.h
HMesh/off_save.o: HMesh/off_save.h CGLA/Vec3d.h CGLA/ArithVec3Float.h
HMesh/off_save.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h CGLA/CGLA-util.h
HMesh/off_save.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
HMesh/off_save.o: CGLA/Vec3usi.h CGLA/Vec3f.h HMesh/Manifold.h
HMesh/off_save.o: HMesh/ConnectivityKernel.h HMesh/ItemVector.h
HMesh/off_save.o: HMesh/ItemID.h HMesh/Iterators.h HMesh/Walker.h
HMesh/off_save.o: HMesh/AttributeVector.h
HMesh/ply_load.o: Geometry/ply_load.h Geometry/TriMesh.h CGLA/Mat4x4f.h
HMesh/ply_load.o: CGLA/ExceptionStandard.h CGLA/CGLA-util.h CGLA/Vec3f.h
HMesh/ply_load.o: CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h CGLA/ArithVec.h
HMesh/ply_load.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
HMesh/ply_load.o: CGLA/Vec3usi.h CGLA/Vec4f.h CGLA/ArithVec4Float.h
HMesh/ply_load.o: CGLA/ArithSqMat4x4Float.h CGLA/ArithSqMatFloat.h
HMesh/ply_load.o: CGLA/ArithMatFloat.h Geometry/IndexedFaceSet.h CGLA/Vec3f.h
HMesh/ply_load.o: CGLA/Vec3i.h Geometry/Material.h Geometry/ply_load.h
HMesh/ply_load.o: HMesh/Manifold.h CGLA/Vec3d.h HMesh/ConnectivityKernel.h
HMesh/ply_load.o: HMesh/ItemVector.h HMesh/ItemID.h HMesh/Iterators.h
HMesh/ply_load.o: HMesh/Walker.h HMesh/AttributeVector.h
HMesh/polygonize.o: Geometry/Implicit.h CGLA/Vec3d.h CGLA/ArithVec3Float.h
HMesh/polygonize.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h CGLA/CGLA-util.h
HMesh/polygonize.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
HMesh/polygonize.o: CGLA/Vec3usi.h CGLA/Vec3f.h Geometry/RGrid.h CGLA/Vec3i.h
HMesh/polygonize.o: Geometry/AncestorGrid.h Geometry/XForm.h
HMesh/polygonize.o: HMesh/polygonize.h Geometry/XForm.h Geometry/RGrid.h
HMesh/polygonize.o: HMesh/Manifold.h HMesh/ConnectivityKernel.h
HMesh/polygonize.o: HMesh/ItemVector.h HMesh/ItemID.h HMesh/Iterators.h
HMesh/polygonize.o: HMesh/Walker.h HMesh/AttributeVector.h HMesh/smooth.h
HMesh/polygonize.o: HMesh/cleanup.h HMesh/triangulate.h
HMesh/quadric_simplify.o: HMesh/quadric_simplify.h HMesh/Manifold.h
HMesh/quadric_simplify.o: CGLA/Vec3d.h CGLA/ArithVec3Float.h
HMesh/quadric_simplify.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h
HMesh/quadric_simplify.o: CGLA/CGLA-util.h CGLA/Vec3i.h CGLA/ArithVec3Int.h
HMesh/quadric_simplify.o: CGLA/ArithVecInt.h CGLA/Vec3usi.h CGLA/Vec3f.h
HMesh/quadric_simplify.o: HMesh/ConnectivityKernel.h HMesh/ItemVector.h
HMesh/quadric_simplify.o: HMesh/ItemID.h HMesh/Iterators.h HMesh/Walker.h
HMesh/quadric_simplify.o: HMesh/AttributeVector.h Geometry/QEM.h
HMesh/quadric_simplify.o: CGLA/Mat3x3d.h CGLA/Vec3d.h
HMesh/quadric_simplify.o: CGLA/ArithSqMat3x3Float.h CGLA/ExceptionStandard.h
HMesh/quadric_simplify.o: CGLA/ArithSqMatFloat.h CGLA/ArithMatFloat.h
HMesh/quadric_simplify.o: HMesh/smooth.h
HMesh/refine_edges.o: HMesh/refine_edges.h HMesh/Manifold.h CGLA/Vec3d.h
HMesh/refine_edges.o: CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h
HMesh/refine_edges.o: CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec3i.h
HMesh/refine_edges.o: CGLA/ArithVec3Int.h CGLA/ArithVecInt.h CGLA/Vec3usi.h
HMesh/refine_edges.o: CGLA/Vec3f.h HMesh/ConnectivityKernel.h
HMesh/refine_edges.o: HMesh/ItemVector.h HMesh/ItemID.h HMesh/Iterators.h
HMesh/refine_edges.o: HMesh/Walker.h HMesh/AttributeVector.h
HMesh/smooth.o: HMesh/smooth.h HMesh/Manifold.h CGLA/Vec3d.h
HMesh/smooth.o: CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h CGLA/ArithVec.h
HMesh/smooth.o: CGLA/CGLA-util.h CGLA/Vec3i.h CGLA/ArithVec3Int.h
HMesh/smooth.o: CGLA/ArithVecInt.h CGLA/Vec3usi.h CGLA/Vec3f.h
HMesh/smooth.o: HMesh/ConnectivityKernel.h HMesh/ItemVector.h HMesh/ItemID.h
HMesh/smooth.o: HMesh/Iterators.h HMesh/Walker.h HMesh/AttributeVector.h
HMesh/smooth.o: CGLA/Mat3x3d.h CGLA/Vec3d.h CGLA/ArithSqMat3x3Float.h
HMesh/smooth.o: CGLA/ExceptionStandard.h CGLA/ArithSqMatFloat.h
HMesh/smooth.o: CGLA/ArithMatFloat.h CGLA/Quatd.h CGLA/ArithQuat.h
HMesh/smooth.o: CGLA/Vec4d.h CGLA/ArithVec4Float.h CGLA/Mat3x3d.h
HMesh/smooth.o: CGLA/Mat4x4d.h CGLA/ArithSqMat4x4Float.h CGLA/Vec4f.h
HMesh/smooth.o: Util/Timer.h
HMesh/subdivision.o: HMesh/dual.h HMesh/subdivision.h CGLA/Vec3d.h
HMesh/subdivision.o: CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h
HMesh/subdivision.o: CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec3i.h
HMesh/subdivision.o: CGLA/ArithVec3Int.h CGLA/ArithVecInt.h CGLA/Vec3usi.h
HMesh/subdivision.o: CGLA/Vec3f.h HMesh/Manifold.h HMesh/ConnectivityKernel.h
HMesh/subdivision.o: HMesh/ItemVector.h HMesh/ItemID.h HMesh/Iterators.h
HMesh/subdivision.o: HMesh/Walker.h HMesh/AttributeVector.h
HMesh/triangulate.o: HMesh/triangulate.h HMesh/Manifold.h CGLA/Vec3d.h
HMesh/triangulate.o: CGLA/ArithVec3Float.h CGLA/ArithVecFloat.h
HMesh/triangulate.o: CGLA/ArithVec.h CGLA/CGLA-util.h CGLA/Vec3i.h
HMesh/triangulate.o: CGLA/ArithVec3Int.h CGLA/ArithVecInt.h CGLA/Vec3usi.h
HMesh/triangulate.o: CGLA/Vec3f.h HMesh/ConnectivityKernel.h
HMesh/triangulate.o: HMesh/ItemVector.h HMesh/ItemID.h HMesh/Iterators.h
HMesh/triangulate.o: HMesh/Walker.h HMesh/AttributeVector.h
HMesh/x3d_load.o: HMesh/x3d_load.h CGLA/Vec3f.h CGLA/ArithVec3Float.h
HMesh/x3d_load.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h CGLA/CGLA-util.h
HMesh/x3d_load.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
HMesh/x3d_load.o: CGLA/Vec3usi.h Util/Timer.h Util/Parse.h CGLA/Vec4f.h
HMesh/x3d_load.o: CGLA/Vec3f.h CGLA/ArithVec4Float.h CGLA/Vec2f.h
HMesh/x3d_load.o: CGLA/ArithVec2Float.h CGLA/Vec2i.h Util/XmlParser.h
HMesh/x3d_load.o: HMesh/Manifold.h CGLA/Vec3d.h HMesh/ConnectivityKernel.h
HMesh/x3d_load.o: HMesh/ItemVector.h HMesh/ItemID.h HMesh/Iterators.h
HMesh/x3d_load.o: HMesh/Walker.h HMesh/AttributeVector.h Geometry/load.h
HMesh/x3d_load.o: Geometry/TriMesh.h
HMesh/x3d_save.o: HMesh/x3d_save.h CGLA/Vec3f.h CGLA/ArithVec3Float.h
HMesh/x3d_save.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h CGLA/CGLA-util.h
HMesh/x3d_save.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
HMesh/x3d_save.o: CGLA/Vec3usi.h HMesh/Manifold.h CGLA/Vec3d.h CGLA/Vec3f.h
HMesh/x3d_save.o: HMesh/ConnectivityKernel.h HMesh/ItemVector.h
HMesh/x3d_save.o: HMesh/ItemID.h HMesh/Iterators.h HMesh/Walker.h
HMesh/x3d_save.o: HMesh/AttributeVector.h
Util/Parse.o: Util/Parse.h CGLA/Vec4f.h CGLA/Vec3f.h CGLA/ArithVec3Float.h
Util/Parse.o: CGLA/ArithVecFloat.h CGLA/ArithVec.h CGLA/CGLA-util.h
Util/Parse.o: CGLA/Vec3i.h CGLA/ArithVec3Int.h CGLA/ArithVecInt.h
Util/Parse.o: CGLA/Vec3usi.h CGLA/ArithVec4Float.h CGLA/Vec3f.h CGLA/Vec2f.h
Util/Parse.o: CGLA/ArithVec2Float.h CGLA/Vec2i.h
Util/XmlParser.o: Util/XmlParser.h Util/string_utils.h
Util/string_utils.o: Util/string_utils.h
