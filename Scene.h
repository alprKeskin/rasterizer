#ifndef _SCENE_H_
#define _SCENE_H_
#include "Vec3.h"
#include "Vec4.h"
#include "Color.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	std::vector<std::vector<Color> > image;
	std::vector<std::vector<double> > depth;
	std::vector<Camera *> cameras;
	std::vector<Vec3 *> vertices;
	std::vector<Color *> colorsOfVertices;
	std::vector<Scaling *> scalings;
	std::vector<Rotation *> rotations;
	std::vector<Translation *> translations;
	std::vector<Mesh *> meshes;

	Scene(const char *xmlPath);

	void assignColorToPixel(int i, int j, Color c);
	void initializeImage(Camera *camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera *camera);
	void convertPPMToPNG(std::string ppmFileName);
	void forwardRenderingPipeline(Camera *camera);
	void lineRasterization(Vec3 v1, Vec3 v2, Camera* camera, Color& c1, Color& c2);
	void triangleRasterization(Camera* camera, Vec3 v1, Vec3 v2, Vec3 v3);
	Color interpolate(double alpha, double beta, double gamma, Vec3 v1, Vec3 v2, Vec3 v3);
	bool liangBarsky(Vec4& firstVec, Vec4& secondVec, Color& color1, Color& color2);

};

#endif
