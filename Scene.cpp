
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"

using namespace tinyxml2;
using namespace std;

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

void Scene::assignColorToPixel(int i, int j, Color c)
{
	this->image[i][j].r = c.r;
	this->image[i][j].g = c.g;
	this->image[i][j].b = c.b;
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;
			vector<double> rowOfDepths;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
				rowOfDepths.push_back(1.01);
			}

			this->image.push_back(rowOfColors);
			this->depth.push_back(rowOfDepths);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				assignColorToPixel(i, j, this->backgroundColor);
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
*/
void Scene::convertPPMToPNG(string ppmFileName)
{
	string command;

	// TODO: Change implementation if necessary.
	command = "./magick convert " + ppmFileName + " " + ppmFileName + ".png";
	system(command.c_str());
}

// ====================================== Our Codes ==========================================

Matrix4 addMatrices(const Matrix4 matrix1, const Matrix4 matrix2) {
	Matrix4 result;
	
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			result.values[i][j] = matrix1.values[i][j] + matrix2.values[i][j];
		}
	}
	return result;
}

Matrix4 multiplyMatrixByScalar(const Matrix4 matrix, const double scaleFactor) {
	Matrix4 result;
	
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			result.values[i][j] = scaleFactor * matrix.values[i][j];
		}
	}
	return result;
}

Color multiplyColorByScalar(Color color, int scalar) {
	return Color(color.r * scalar, color.g * scalar, color.b * scalar);
}

Color subtractColors(Color c1, Color c2) {
	return Color(c1.r - c2.r, c1.g - c2.g, c1.b - c2.b);
}

Color addColors(Color c1, Color c2) {
	return Color(c1.r + c2.r, c1.g + c2.g, c1.b + c2.b);
}

Color divideColorByScalar(Color color, int division) {
	return Color(color.r / division, color.g / division, color.b / division);
}

Matrix4 getTranslationMatrix(const Translation* translation) {
	Matrix4 translationMatrix = getIdentityMatrix();
	translationMatrix.values[0][3] = translation -> tx;
	translationMatrix.values[1][3] = translation -> ty;
	translationMatrix.values[2][3] = translation -> tz;		
	return translationMatrix;
}

Matrix4 getScaleMatrix(const Scaling* scaling) {
	Matrix4 scaleMatrix = getIdentityMatrix();
	scaleMatrix.values[0][0] = scaling -> sx;
	scaleMatrix.values[1][1] = scaling -> sy;
	scaleMatrix.values[2][2] = scaling -> sz;
	scaleMatrix.values[3][3] = 1;
	return scaleMatrix;
}

Matrix4 getRotationMatrix(const Rotation* rotation) {
	Matrix4 matrix = getIdentityMatrix();
    Vec3 rotationVector = normalizeVec3(Vec3(rotation->ux, rotation->uy, rotation->uz));
    double ux = rotationVector.x;
    double uy = rotationVector.y;
    double uz = rotationVector.z;
    double radian = (rotation->angle) * M_PI / 180.0;
    double cosA = cos(radian);
    double sinA = sin(radian);

    matrix.values[0][0] = cosA + ux * ux * (1 - cosA);
    matrix.values[0][1] = ux * uy * (1 - cosA) - uz * sinA;
    matrix.values[0][2] = ux * uz * (1 - cosA) + uy * sinA;
    matrix.values[0][3] = 0;
    matrix.values[1][0] = uy * ux * (1 - cosA) + uz * sinA;
    matrix.values[1][1] = cosA + uy * uy * (1 - cosA);
    matrix.values[1][2] = uy * uz * (1 - cosA) - ux * sinA;
    matrix.values[1][3] = 0;
    matrix.values[2][0] = uz * ux * (1 - cosA) - uy * sinA;
    matrix.values[2][1] = uz * uy * (1 - cosA) + ux * sinA;
    matrix.values[2][2] = cosA + uz * uz * (1 - cosA);
    matrix.values[2][3] = 0;
    matrix.values[3][0] = 0;
    matrix.values[3][1] = 0;
    matrix.values[3][2] = 0;
    matrix.values[3][3] = 1;
	return matrix;
}

Matrix4 getModelingMatrixForMesh(Mesh* mesh, vector<Translation*> translations, vector<Scaling*>scalings, vector<Rotation*> rotations) {
	Matrix4 modelingMatrix = getIdentityMatrix();
	
	for (int i = 0; i < mesh -> numberOfTransformations; i++) {
		if (mesh -> transformationTypes[i] == 't') {	
			modelingMatrix = multiplyMatrixWithMatrix(getTranslationMatrix(translations[mesh -> transformationIds[i] - 1]), modelingMatrix);
		}
		else if (mesh -> transformationTypes[i] == 's') {
			modelingMatrix = multiplyMatrixWithMatrix(getScaleMatrix(scalings[mesh -> transformationIds[i] - 1]),modelingMatrix);
		}
		else if (mesh -> transformationTypes[i] == 'r') {
			modelingMatrix = multiplyMatrixWithMatrix(getRotationMatrix(rotations[mesh -> transformationIds[i] - 1]), modelingMatrix);
		}
	}
	return modelingMatrix;
}

Vec4 convertVec3ToVec4(const Vec3& point) {
	return Vec4(point.x, point.y, point.z, 1, point.colorId);
}

Matrix4 cameraTransformations(const Camera& camera){
	Matrix4 cameraMatrix = getIdentityMatrix();

    cameraMatrix.values[0][0] = camera.u.x;
    cameraMatrix.values[0][1] = camera.u.y;
    cameraMatrix.values[0][2] = camera.u.z;
    cameraMatrix.values[0][3] = -1 * ((camera.u.x) * (camera.position.x) + (camera.u.y) * (camera.position.y) + (camera.u.z) * (camera.position.z));

    cameraMatrix.values[1][0] = camera.v.x;
    cameraMatrix.values[1][1] = camera.v.y;
    cameraMatrix.values[1][2] = camera.v.z;
    cameraMatrix.values[1][3] = -1 * ((camera.v.x) * (camera.position.x) + (camera.v.y) * (camera.position.y) + (camera.v.z) * (camera.position.z));
    
	cameraMatrix.values[2][0] = camera.w.x;
    cameraMatrix.values[2][1] = camera.w.y;
    cameraMatrix.values[2][2] = camera.w.z;
    cameraMatrix.values[2][3] = -1 * ((camera.w.x) * (camera.position.x) + (camera.w.y) * (camera.position.y) + (camera.w.z) * (camera.position.z));

    cameraMatrix.values[3][0] = 0.0;
    cameraMatrix.values[3][1] = 0.0;
    cameraMatrix.values[3][2] = 0.0;
    cameraMatrix.values[3][3] = 1.0;

	return cameraMatrix;
}

Matrix4 getOrthographicProjectionMatrix(Camera camera) {
	Matrix4 orthographicProjectionMatrix = getIdentityMatrix();
	
	orthographicProjectionMatrix.values[0][0] = 2 / (camera.right - camera.left);
	orthographicProjectionMatrix.values[1][1] = 2 / (camera.top - camera.bottom);
	orthographicProjectionMatrix.values[2][2] = -2 / (camera.far - camera.near);
	orthographicProjectionMatrix.values[0][3] = -1 * ((camera.right + camera.left) / (camera.right - camera.left));
	orthographicProjectionMatrix.values[1][3] = -1 * ((camera.top + camera.bottom) / (camera.top - camera.bottom));
	orthographicProjectionMatrix.values[2][3] = -1 * ((camera.far + camera.near) / (camera.far - camera.near));
	
	return orthographicProjectionMatrix;
}

Matrix4 getPerspectiveProjectionMatrix(Camera camera) {
	Matrix4 perspectiveProjectionMatrix = getIdentityMatrix();
	
	perspectiveProjectionMatrix.values[0][0] = 2 * camera.near / (camera.right - camera.left);
	perspectiveProjectionMatrix.values[0][2] = (camera.right + camera.left) / (camera.right - camera.left);
	perspectiveProjectionMatrix.values[1][1] = 2 * camera.near / (camera.top - camera.bottom);
	perspectiveProjectionMatrix.values[1][2] = (camera.top + camera.bottom) / (camera.top - camera.bottom);
	perspectiveProjectionMatrix.values[2][2] = -1 * (camera.far + camera.near) / (camera.far - camera.near);
	perspectiveProjectionMatrix.values[3][2] = -1;
	perspectiveProjectionMatrix.values[2][3] = -2 * camera.far * camera.near / (camera.far - camera.near);
	perspectiveProjectionMatrix.values[3][3] = 0;
	
	return perspectiveProjectionMatrix;
}

Vec4 perspectiveDivision(Vec4 v){ 
    return Vec4(v.x/v.t, v.y/v.t, v.z/v.t, 1, v.colorId); 
}

bool isVisible(double den, double num, double& te, double& tl){
    if (den > 0){
        double t = num/den;
        if (t > tl){
            return false;
        }
        if (t > te){
            te = t;
        }
    } 
    else if (den < 0){
        double t = num/den;
        if (t < te){
            return false;
        }
        if (t < tl){
            tl = t;
        }
    } else if (num > 0){
        return false;
    }
    return true;
}

bool Scene::liangBarsky(Vec4& firstVec, Vec4& secondVec, Color& color1, Color& color2){
    double bound = 1.0;
    double te = 0.0;
    double tl = 1.0;
    double dx = secondVec.x-firstVec.x;
    double dy = secondVec.y-firstVec.y;
    double dz = secondVec.z-firstVec.z;
	color1 = *(this->colorsOfVertices[firstVec.colorId - 1]);
	color2 = *(this->colorsOfVertices[secondVec.colorId - 1]);
	Color dColor = divideColorByScalar(subtractColors(color2, color1), 1);
    bool visible = false;
    if (isVisible(dx, -bound-firstVec.x, te, tl)) { //left
        if (isVisible(-dx, firstVec.x-bound, te, tl)) { //right
            if (isVisible(dy, -bound-firstVec.y, te, tl)) { //bottom
                if (isVisible(-dy, firstVec.y-bound, te, tl)) { //top
                    if (isVisible(dz, -bound-firstVec.z, te, tl)) { //front
						if (isVisible(-dz, firstVec.z-bound, te, tl)) { //back
							visible = true;
							if (tl < 1) {
								secondVec.x = firstVec.x+dx*tl;
								secondVec.y = firstVec.y+dy*tl;
								secondVec.z = firstVec.z+dz*tl;
								color2.r = color1.r + dColor.r*tl;
								color2.g = color1.g + dColor.g*tl;
								color2.b = color1.b + dColor.b*tl;
							}
							if (te > 0) {
								firstVec.x = firstVec.x+dx*te;
								firstVec.y = firstVec.y+dy*te;
								firstVec.z = firstVec.z+dz*te;
								color1.r = color1.r + dColor.r*te;
								color1.g = color1.g + dColor.g*te;
								color1.b = color1.b + dColor.b*te;
							}
						}
                    }
                }
            }
        }
    }
    return visible;
}

Vec3 viewportTransform(Camera* cam, Vec4 vector){
    double x = (cam->horRes*(vector.x+1)-1)/2.0;
    double y = (cam->verRes*(vector.y+1)-1)/2.0;
    double z = (vector.z+1)/2.0;
    return Vec3(x,y,z,vector.colorId);
}

void Scene::lineRasterization(Vec3 v1, Vec3 v2, Camera* camera, Color& color1, Color& color2) {
    // Determine the differences in the x, y coordinates
    double deltaX = abs(round(v2.x) - round(v1.x));
    double deltaY = abs(round(v2.y) - round(v1.y));

	// initialize start and end points
	int xStart, yStart, xEnd, yEnd;
	Color startColor, endColor;

    // Retrieve the colors associated with the vertices
    Color colorOfFirstVertex = color1;
    Color colorOfSecondVertex = color2;
	
	if ( ((deltaX > deltaY) && (v2.x >= v1.x)) || ((deltaY >= deltaX) && (v2.y >= v1.y))) {
		xStart = v1.x;
		xEnd = v2.x;
		yStart = v1.y;
		yEnd = v2.y;
		startColor = colorOfFirstVertex;
		endColor = colorOfSecondVertex;	
	}
	else if (((deltaX > deltaY) && (v1.x > v2.x)) || ((deltaY >= deltaX) && (v1.y > v2.y))) {
		xStart = v2.x;
		xEnd = v1.x;
		yStart = v2.y;
		yEnd = v1.y;
		startColor = colorOfSecondVertex;
		endColor = colorOfFirstVertex;	
	}
	else{
		cout << "Function should not be here";
		exit (1);
	}
	int dx = xEnd - xStart;
	int dy = yEnd - yStart;
	
	if (deltaX > deltaY) {
		Color interpolatedColor = startColor; 
		int yInitial = 1;
		if (dy < 0) {
			yInitial = -1;
			dy = -dy;
		}

		Color dColor = divideColorByScalar(subtractColors(endColor, startColor), dx);
		int dMultiplied = 2*dy-dx;
		int y = yStart;

		for (int x = xStart; x < xEnd; x++) {
			double zValue = (double) (x - xStart) / (xEnd - xStart) * (v2.z - v1.z) + v1.z;
			
			if (x < camera->horRes && y < camera->verRes && x >= 0 && y >= 0 && zValue < this->depth[x][y]){
				this->image[x][y] = interpolatedColor;
				this->depth[x][y] = zValue;
			}

			if (dMultiplied > 0) {
				y += yInitial;
				dMultiplied += 2*(dy-dx);
			} else {
				dMultiplied += 2*dy;
			}
			interpolatedColor = addColors(interpolatedColor, dColor);
		}
	}

	else if (deltaY >= deltaX) {
		Color interpolatedColor = startColor;
		
		int xInitial = 1;
		if (dx < 0) {
			xInitial = -1;
			dx = -dx;
		}

		Color dColor = divideColorByScalar(subtractColors(endColor, startColor), dy);
		int dMultiplied = 2*dx-dy;
		int x = xStart;

		for (int y = yStart; y < yEnd; y++) {
			// mark the point with color c.
			if (x < camera->horRes && y < camera->verRes && x >= 0 && y >= 0){
				this->image[x][y] = interpolatedColor;
			}

			if (dMultiplied > 0) {
				x += xInitial;
				dMultiplied += 2*(dx-dy);
			} else {
				dMultiplied += 2*dx;
			}
			interpolatedColor = addColors(interpolatedColor, dColor);
		}
	}
}

double edgeFunc(const Vec3& a, const Vec3& b, const Vec3& c) {
    return (c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x);
}

void Scene::triangleRasterization(Camera *camera, Vec3 v1, Vec3 v2, Vec3 v3) {

    int minX = round(min(min(v1.x, v2.x), v3.x)), maxX = round(max(max(v1.x, v2.x), v3.x));
    int minY = round(min(min(v1.y, v2.y), v3.y)), maxY = round(max(max(v1.y, v2.y), v3.y));

    for (int y = minY; y <= maxY; ++y) {
        for (int x = minX; x <= maxX; ++x) {
            Vec3 p(x, y, 0);
            double alpha = edgeFunc(v2, v3, p) / edgeFunc(v2, v3, v1);
            double beta  = edgeFunc(v3, v1, p) / edgeFunc(v3, v1, v2);
            double gamma = edgeFunc(v1, v2, p) / edgeFunc(v1, v2, v3);

            if (alpha >= 0 && beta >= 0 && gamma >= 0) {
                p.z = alpha * v1.z + beta * v2.z + gamma * v3.z;
                if (x >= 0 && y >= 0 && x < camera->horRes && y < camera->verRes && p.z < this->depth[x][y]) {
                    this->image[x][y] = interpolate(alpha, beta, gamma, v1, v2, v3);
                    this->depth[x][y] = p.z;
                }
            }
        }
    }
}

Color Scene::interpolate(double alpha, double beta, double gamma, Vec3 v1, Vec3 v2, Vec3 v3) {
    return Color(
        alpha * this->colorsOfVertices[v1.colorId - 1]->r + beta * this->colorsOfVertices[v2.colorId - 1]->r + gamma * this->colorsOfVertices[v3.colorId - 1]->r,
        alpha * this->colorsOfVertices[v1.colorId - 1]->g + beta * this->colorsOfVertices[v2.colorId - 1]->g + gamma * this->colorsOfVertices[v3.colorId - 1]->g,
        alpha * this->colorsOfVertices[v1.colorId - 1]->b + beta * this->colorsOfVertices[v2.colorId - 1]->b + gamma * this->colorsOfVertices[v3.colorId - 1]->b
    );
}

// Transformations, clipping, culling, rasterization are done here.
void Scene::forwardRenderingPipeline(Camera *camera) {
	vector<Vec3*> vertices = this -> vertices; 
	vector<Mesh*> meshes = this -> meshes;
	vector<Translation*> translations = this -> translations;
	vector<Scaling*> scalings = this -> scalings;
	vector<Rotation*> rotations = this -> rotations;

	Matrix4 modellingMatrix;

	Matrix4 matrixCam = cameraTransformations(*camera);
	for (int i = 0; i < meshes.size(); i++) {
		Mesh* mesh = meshes[i];

		Matrix4 matrixModel = getModelingMatrixForMesh(mesh, translations, scalings, rotations);
		modellingMatrix = multiplyMatrixWithMatrix(matrixCam, matrixModel);
		
		if(camera -> projectionType == ORTOGRAPHIC_PROJECTION) {
			modellingMatrix = multiplyMatrixWithMatrix(getOrthographicProjectionMatrix(*camera), modellingMatrix);
		}
		else {
			modellingMatrix = multiplyMatrixWithMatrix(getPerspectiveProjectionMatrix(*camera),modellingMatrix);
		}
		
		for (int k = 0; k < mesh -> numberOfTriangles; k++) {
			Triangle triangle = mesh -> triangles[k];
			
			Vec3 vertex0 = *(this -> vertices[triangle.vertexIds[0] - 1]);
			Vec3 vertex1 = *(this -> vertices[triangle.vertexIds[1] - 1]);
			Vec3 vertex2 = *(this -> vertices[triangle.vertexIds[2] - 1]);

			Vec4 transformedVertex0 = multiplyMatrixWithVec4(modellingMatrix, convertVec3ToVec4(vertex0));
			Vec4 transformedVertex1 = multiplyMatrixWithVec4(modellingMatrix, convertVec3ToVec4(vertex1));
			Vec4 transformedVertex2 = multiplyMatrixWithVec4(modellingMatrix, convertVec3ToVec4(vertex2));
			
			transformedVertex0 = perspectiveDivision(transformedVertex0);
			transformedVertex1 = perspectiveDivision(transformedVertex1);
			transformedVertex2 = perspectiveDivision(transformedVertex2);

			// Culling
			if (this -> cullingEnabled) {
				Vec3 vector1 = Vec3(transformedVertex1.x - transformedVertex0.x, transformedVertex1.y - transformedVertex0.y, transformedVertex1.z - transformedVertex0.z);
				Vec3 vector2 = Vec3(transformedVertex2.x - transformedVertex0.x, transformedVertex2.y - transformedVertex0.y, transformedVertex2.z - transformedVertex0.z);
				Vec3 normal = normalizeVec3(crossProductVec3(vector2, vector1));
				Vec3 sumOfVertices = Vec3(transformedVertex0.x + transformedVertex1.x + transformedVertex2.x, transformedVertex0.y + transformedVertex1.y + transformedVertex2.y, transformedVertex0.z + transformedVertex1.z + transformedVertex2.z);
				// center of triangle
				Vec3 cameraVector = Vec3(sumOfVertices.x / 3, sumOfVertices.y / 3, sumOfVertices.z / 3);
				if (dotProductVec3(normal, cameraVector) >= 0) {
					// cull it
					continue;
				}
			}
			// Clipping
			if (mesh -> type == WIREFRAME_MESH) {
				// Line Rasterization
				Vec4 edge1To1 = transformedVertex0;
                Vec4 edge1To2 = transformedVertex1;
                Vec4 edge2To1 = transformedVertex1;
                Vec4 edge2To2 = transformedVertex2;
                Vec4 edge3To1 = transformedVertex0;
                Vec4 edge3To2 = transformedVertex2;

				Color color1, color2, color3, color4, color5, color6;

                bool edge1Visible = liangBarsky(edge1To1, edge1To2, color1, color2);
                bool edge2Visible = liangBarsky(edge2To1, edge2To2, color3, color4);
                bool edge3Visible = liangBarsky(edge3To1, edge3To2, color5, color6);
				/*
				cout << "color1: " << color1 . r << " " << color1 . g << " " << color1 . b << endl;
				cout << "color2: " << color2 . r << " " << color2 . g << " " << color2 . b << endl;
				cout << "color3: " << color3 . r << " " << color3 . g << " " << color3 . b << endl;
				cout << "color4: " << color4 . r << " " << color4 . g << " " << color4 . b << endl;
				cout << "color5: " << color5 . r << " " << color5 . g << " " << color5 . b << endl;
				cout << "color6: " << color6 . r << " " << color6 . g << " " << color6 . b << endl;
                */
				if (edge1Visible) {
                    Vec3 edge1To1Viewported = viewportTransform(camera, edge1To1);
                    Vec3 edge1To2Viewported = viewportTransform(camera, edge1To2);
                    lineRasterization(edge1To1Viewported, edge1To2Viewported, camera, color1, color2);
                }
                if (edge2Visible) {
                    Vec3 edge2To1Viewported = viewportTransform(camera, edge2To1);
                    Vec3 edge2To2Viewported = viewportTransform(camera, edge2To2);
                    lineRasterization(edge2To1Viewported, edge2To2Viewported, camera, color3, color4);
                }
                if (edge3Visible) {
                    Vec3 edge3To1Viewported = viewportTransform(camera, edge3To1);
                    Vec3 edge3To2Viewported = viewportTransform(camera, edge3To2);
                    lineRasterization(edge3To1Viewported, edge3To2Viewported, camera, color5, color6);
                }
			}
			else {
				// Triangle Rasterization
				Vec3 vecView1 = viewportTransform(camera, transformedVertex0);
				Vec3 vecView2 = viewportTransform(camera, transformedVertex1);
				Vec3 vecView3 = viewportTransform(camera, transformedVertex2);
				triangleRasterization(camera, vecView1, vecView2, vecView3);
			}
		}

	}
}
