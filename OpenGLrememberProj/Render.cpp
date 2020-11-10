#include "Render.h"
#include "Vec.h"

#include <sstream>
#include <iostream>
#include <Array>
#include <array>
#include <vector>
using std::vector;
using std::array;


#include <windows.h>
#include <GL\GL.h>
#include <GL\GLU.h>

#include "MyOGL.h"

#include "Camera.h"
#include "Light.h"
#include "Primitives.h"

#include "GUItextRectangle.h"

bool textureMode = true;
bool lightMode = true;
bool alphaMode = false;

//класс для настройки камеры
class CustomCamera : public Camera
{
public:
	//дистанция камеры
	double camDist;
	//углы поворота камеры
	double fi1, fi2;

	
	//значния масеры по умолчанию
	CustomCamera()
	{
		camDist = 15;
		fi1 = 1;
		fi2 = 1;
	}

	
	//считает позицию камеры, исходя из углов поворота, вызывается движком
	void SetUpCamera()
	{
		//отвечает за поворот камеры мышкой
		lookPoint.setCoords(0, 0, 0);

		pos.setCoords(camDist*cos(fi2)*cos(fi1),
			camDist*cos(fi2)*sin(fi1),
			camDist*sin(fi2));

		if (cos(fi2) <= 0)
			normal.setCoords(0, 0, -1);
		else
			normal.setCoords(0, 0, 1);

		LookAt();
	}

	void CustomCamera::LookAt()
	{
		//функция настройки камеры
		gluLookAt(pos.X(), pos.Y(), pos.Z(), lookPoint.X(), lookPoint.Y(), lookPoint.Z(), normal.X(), normal.Y(), normal.Z());
	}



}  camera;   //создаем объект камеры


//Класс для настройки света
class CustomLight : public Light
{
public:
	CustomLight()
	{
		//начальная позиция света
		pos = Vector3(1, 1, 3);
	}

	
	//рисует сферу и линии под источником света, вызывается движком
	void  DrawLightGhismo()
	{
		glDisable(GL_LIGHTING);

		glColor3d(0.9, 0.8, 0);
		Sphere s;
		s.pos = pos;
		s.scale = s.scale*0.08;
		s.Show();
		
		if (OpenGL::isKeyPressed('G'))
		{
			glColor3d(0, 0, 0);
			//линия от источника света до окружности
			glBegin(GL_LINES);
			glVertex3d(pos.X(), pos.Y(), pos.Z());
			glVertex3d(pos.X(), pos.Y(), 0);
			glEnd();

			//рисуем окруность
			Circle c;
			c.pos.setCoords(pos.X(), pos.Y(), 0);
			c.scale = c.scale*1.5;
			c.Show();
		}

	}

	void SetUpLight()
	{
		GLfloat amb[] = { 0.2, 0.2, 0.2, 0 };
		GLfloat dif[] = { 1.0, 1.0, 1.0, 0 };
		GLfloat spec[] = { .7, .7, .7, 0 };
		GLfloat position[] = { pos.X(), pos.Y(), pos.Z(), 1. };

		// параметры источника света
		glLightfv(GL_LIGHT0, GL_POSITION, position);
		// характеристики излучаемого света
		// фоновое освещение (рассеянный свет)
		glLightfv(GL_LIGHT0, GL_AMBIENT, amb);
		// диффузная составляющая света
		glLightfv(GL_LIGHT0, GL_DIFFUSE, dif);
		// зеркально отражаемая составляющая света
		glLightfv(GL_LIGHT0, GL_SPECULAR, spec);

		glEnable(GL_LIGHT0);
	}


} light;  //создаем источник света

std::vector<std::vector<Vec>> points({
	{
		{0, 0, 3},
		{0, -1, 2},
		{0, -2, 2},
		{0, -3, 1}
	},
	{
		{1, 0, 1},
		{1, -1, 2},
		{1, -2, -2},
		{1, -3, 1}
	},
	{
		{2, 0, 1},
		{2, -1, 2},
		{2, -2, 2},
		{2, -3, 1}
	},
	{
		{3, 0, 1},
		{3, -1, 2},
		{3, -2, 2},
		{3, -3, 1}
	}
	});

void find_and_change_point(LPPOINT poi, double dy) {
	GLint    viewport[4];    // параметры viewport-a.
	GLdouble projection[16]; // матрица проекции.
	GLdouble modelview[16];  // видовая матрица.

	glGetIntegerv(GL_VIEWPORT, viewport);           // узнаём параметры viewport-a.
	glGetDoublev(GL_PROJECTION_MATRIX, projection); // узнаём матрицу проекции.
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);   // узнаём видовую матрицу.
	Vec p = { 0, 0, 0 };
	double delta = 10;
	for (auto& v : points) {
		for (auto& elem : v) {
			gluProject(elem.x, elem.y, elem.z, modelview, projection, viewport, &p.x, &p.y, &p.z);
			if (p.x > poi->x - delta && p.x < poi->x + delta &&
				p.y > poi->y - delta && p.y < poi->y + delta) {
				elem.z += 0.02 * dy;
			}
		}
	}
};

//старые координаты мыши
int mouseX = 0, mouseY = 0;

void mouseEvent(OpenGL *ogl, int mX, int mY)
{
	int dx = mouseX - mX;
	int dy = mouseY - mY;
	mouseX = mX;
	mouseY = mY;

	//меняем углы камеры при нажатой левой кнопке мыши
	if (OpenGL::isKeyPressed(VK_RBUTTON))
	{
		camera.fi1 += 0.01*dx;
		camera.fi2 += -0.01*dy;
	}

	
	//двигаем свет по плоскости, в точку где мышь
	if (OpenGL::isKeyPressed('G') && !OpenGL::isKeyPressed(VK_LBUTTON))
	{
		LPPOINT POINT = new tagPOINT();
		GetCursorPos(POINT);
		ScreenToClient(ogl->getHwnd(), POINT);
		POINT->y = ogl->getHeight() - POINT->y;

		Ray r = camera.getLookRay(POINT->x, POINT->y);

		double z = light.pos.Z();

		double k = 0, x = 0, y = 0;
		if (r.direction.Z() == 0)
			k = 0;
		else
			k = (z - r.origin.Z()) / r.direction.Z();

		x = k*r.direction.X() + r.origin.X();
		y = k*r.direction.Y() + r.origin.Y();

		light.pos = Vector3(x, y, z);
	}

	if (OpenGL::isKeyPressed('G') && OpenGL::isKeyPressed(VK_LBUTTON))
	{
		light.pos = light.pos + Vector3(0, 0, 0.02*dy);
	}

	if (OpenGL::isKeyPressed(VK_LBUTTON))
	{
		LPPOINT POINT = new tagPOINT();
		GetCursorPos(POINT);
		ScreenToClient(ogl->getHwnd(), POINT);
		POINT->y = ogl->getHeight() - POINT->y;
		find_and_change_point(POINT, dy);
		delete POINT;
	}

	if (OpenGL::isKeyPressed('G') && OpenGL::isKeyPressed(VK_LBUTTON))
	{
		light.pos = light.pos + Vector3(0, 0, 0.02 * dy);
	}	
}


void mouseWheelEvent(OpenGL *ogl, int delta)
{

	if (delta < 0 && camera.camDist <= 1)
		return;
	if (delta > 0 && camera.camDist >= 100)
		return;

	camera.camDist += 0.01*delta;

}

void keyDownEvent(OpenGL *ogl, int key)
{
	if (key == 'L')
	{
		lightMode = !lightMode;
	}

	if (key == 'T')
	{
		textureMode = !textureMode;
	}

	if (key == 'A')
	{
		alphaMode = !alphaMode;
	}

	if (key == 'R')
	{
		camera.fi1 = 1;
		camera.fi2 = 1;
		camera.camDist = 15;

		light.pos = Vector3(1, 1, 3);
	}

	if (key == 'F')
	{
		light.pos = camera.pos;
	}
}

void keyUpEvent(OpenGL *ogl, int key)
{
	
}



GLuint texId;

//выполняется перед первым рендером
void initRender(OpenGL *ogl)
{
	//настройка текстур

	//4 байта на хранение пикселя
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);

	//настройка режима наложения текстур
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	//включаем текстуры
	glEnable(GL_TEXTURE_2D);
	

	//массив трехбайтных элементов  (R G B)
	RGBTRIPLE *texarray;

	//массив символов, (высота*ширина*4      4, потомучто   выше, мы указали использовать по 4 байта на пиксель текстуры - R G B A)
	char *texCharArray;
	int texW, texH;
	OpenGL::LoadBMP("texture.bmp", &texW, &texH, &texarray);
	OpenGL::RGBtoChar(texarray, texW, texH, &texCharArray);

	
	//генерируем ИД для текстуры
	glGenTextures(1, &texId);
	//биндим айдишник, все что будет происходить с текстурой, будте происходить по этому ИД
	glBindTexture(GL_TEXTURE_2D, texId);

	//загружаем текстуру в видеопямять, в оперативке нам больше  она не нужна
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texW, texH, 0, GL_RGBA, GL_UNSIGNED_BYTE, texCharArray);

	//отчистка памяти
	free(texCharArray);
	delete[] texarray;

	//наводим шмон
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);


	//камеру и свет привязываем к "движку"
	ogl->mainCamera = &camera;
	ogl->mainLight = &light;

	// нормализация нормалей : их длины будет равна 1
	glEnable(GL_NORMALIZE);

	// устранение ступенчатости для линий
	glEnable(GL_LINE_SMOOTH); 


	//   задать параметры освещения
	//  параметр GL_LIGHT_MODEL_TWO_SIDE - 
	//                0 -  лицевые и изнаночные рисуются одинаково(по умолчанию), 
	//                1 - лицевые и изнаночные обрабатываются разными режимами       
	//                соответственно лицевым и изнаночным свойствам материалов.    
	//  параметр GL_LIGHT_MODEL_AMBIENT - задать фоновое освещение, 
	//                не зависящее от сточников
	// по умолчанию (0.2, 0.2, 0.2, 1.0)

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 0);

	camera.fi1 = -1.3;
	camera.fi2 = 0.8;
}
Vec TexPoint(const Vec& tmp) {
	
	return {tmp.x / 3, tmp.y / 3, 0};
}

unsigned long long factorial(int n) {
	unsigned long long factor = 1;
	for (int i = 1; i <= n; ++i) {
		factor *= i;
	}
	return factor;
}

double Bernuly(double u, double n, int index) {
	return (factorial(n) / (factorial(index) * factorial(n - index))) * pow(u, index) * pow(1 - u, n - index);
}

void BezeSurfacePoint(double u, double v, Vec& a) {
	Vec new_v;
	int n = 3, m = 3;
	for (size_t i = 0; i < points.size(); ++i) {
		for (size_t j = 0; j < points[i].size(); ++j) {
			auto fact1 = Bernuly(u, n, i);
			auto fact2 = Bernuly(v, m, j);
			new_v += points[i][j] * fact1 * fact2;
		}
	}
	a = new_v;
}

double* Normal(Vec first, Vec second, Vec third, bool reverse = false) {
	Vec a = { second.x - first.x, second.y - first.y, second.z - first.z };
	Vec b = { third.x - first.x, third.y - first.y, third.z - first.z };
	Vec norm = { a.y * b.z - b.y * a.z, -a.x * b.z + b.x * a.z, a.x * b.y - b.x * a.y };
	double len = sqrt(pow(norm.x, 2) + pow(norm.y, 2) + pow(norm.z, 2));
	if (reverse) {
		Vec tmp = { norm.x / len, norm.y / len, norm.z / len };
		tmp *= -1;
		return (double*)&Vec(tmp);
	}
	return (double*)&Vec(norm.x / len, norm.y / len, norm.z / len);
}

//Vec Normal(const double(&A)[3], const double(&B)[3], const double(&C)[3]) {
//	double A1[3], C1[3];
//	Vec N;
//	for (int i = 0; i < 3; i++) {
//		A1[i] = B[i] - A[i];
//		C1[i] = B[i] - C[i];
//	}
//	N.x = (A1[1] * C1[2]) - (C1[1] * A1[2]);
//	N.y = (C1[0] * A1[2]) - (A1[0] * C1[2]);
//	N.z = (A1[0] * C1[1]) - (C1[0] * A1[1]);
//	double len = sqrt(pow(N.x, 2) + pow(N.y, 2) + pow(N.z, 2));
//		N / len;
//	return N;
//}

Vec F_Beze(Vec p1, Vec p2, Vec p3, Vec p4, double t) {
	return p1 * (1 - t) * (1 - t) * (1 - t) + p2 * 3 * t * (1 - t) * (1 - t) + p3 * t * t * (1 - t) + p4 * t * t * t;
}

Vec F_Ermit(Vec p1, Vec p4, Vec r1, Vec r4, double t) {

	return p1 * (2 * t * t * t - 3 * t * t + 1) + p4 * (-2 * t * t * t + 3 * t * t) + r1 * (t * t * t - 2 * t * t + t) + r4 * (t * t * t - t * t);
}

void Beze(Vec P1, Vec P2, Vec P3, Vec P4, double t_max = 1) {
	Vec P;

	glLineWidth(3);
	glColor3d(0.5, 0, 0.5);

	glDisable(GL_LIGHTING);
	glBegin(GL_LINE_STRIP);
	glVertex3dv(P1.ToArray());
	glVertex3dv(P2.ToArray());
	glVertex3dv(P3.ToArray());
	glVertex3dv(P4.ToArray());
	glEnd();

	glBegin(GL_LINE_STRIP);
	glLineWidth(5);
	glColor3d(0.5, 0.5, 0);
	for (double t = 0; t < t_max; t += 0.01)
	{
		P = F_Beze(P1, P2, P3, P4, t);
		glVertex3dv(P.ToArray());
	}
	glEnd();
	glEnable(GL_LIGHTING);
}

void Ermit(Vec P1, Vec P2, Vec P3, Vec P4, double t_max) {
	Vec P;

	glLineWidth(3);
	glColor3d(0, 0.5, 0.5);

	glBegin(GL_LINES);
	glVertex3dv(P1.ToArray());
	glVertex3dv(P2.ToArray());
	glVertex3dv(P3.ToArray());
	glVertex3dv(P4.ToArray());
	glEnd();

	glBegin(GL_LINE_STRIP);
	glColor3d(0, 0.1, 0.9);
	glLineWidth(5);
	Vec R1;
	Vec R4;
	R1 = P2 - P1;
	R4 = P4 - P3;

	for (double t = 0; t <= t_max; t += 0.01)
	{
		P = F_Ermit(P1, P4, R1, R4, t);
		glVertex3dv(P.ToArray());
	}
	glEnd();
}

void _30(double delta_time, double t_max) {
	Vec P1 = { 0, 0, 0 };
	Vec P2 = { -4, 6, 7 };
	Vec P3 = { 4, 8, 2 };
	Vec P4 = { 10, 10, 0 };

	Beze(P1, P2, P3, P4, t_max);

	Vec A1 = { 3,3,3 };
	Vec A2 = { 4,6,-7 };
	Vec A3 = { 7, -1, 2 };
	Vec A4 = { 0,0,0 };

	Beze(A1, A2, A3, A4, t_max);

	Vec B1 = { -3, 0, 0 };
	Vec B2 = { 500, -300, 7 };
	Vec B3 = { 9, -10, 2 };
	Vec B4 = { 10, -15, 0 };

	Ermit(B1, B2, B3, B4, t_max);

	Vec C1 = { 0, 0, 0 };
	Vec C2 = { -1, -5, 7 };
	Vec C3 = { -3, -8, -2 };
	Vec C4 = { -10, -8, 0 };

	Ermit(C1, C2, C3, C4, t_max);
}

void Square(const Vec& O) {
	Vec bot[4];
	Vec top[4];
	double h = 0.5;

	bot[0].x = O.x + h;//x
	bot[1].x = O.x - h;
	bot[2].x = O.x - h;
	bot[3].x = O.x + h;
					 
	bot[0].y = O.y + h;//y
	bot[1].y = O.y + h;
	bot[2].y = O.y - h;
	bot[3].y = O.y - h;

	for (int i = 0; i < 4; i++) {//z
		bot[i].z = O.z - h;
	}

	for (int i = 0; i < 4; i++) {
		top[i].x = bot[i].x;
		top[i].y = bot[i].y;
		top[i].z = O.z + h;//z
	}

	glBegin(GL_QUAD_STRIP);
	glColor3d(1, 0, 0);
	for (int i = 0; i < 3; i++) {
		glNormal3dv(Normal(top[i], bot[i], top[i + 1]));
		glVertex3dv(top[i].ToArray());
		glVertex3dv(bot[i].ToArray());
		glVertex3dv(top[i + 1].ToArray());
		glVertex3dv(bot[i + 1].ToArray());
	}
	glNormal3dv(Normal(top[3], bot[3], top[0]));
	glVertex3dv(top[3].ToArray());
	glVertex3dv(bot[3].ToArray());
	glVertex3dv(top[0].ToArray());
	glVertex3dv(bot[0].ToArray());
	glEnd();
	glBegin(GL_QUADS);
	glColor3d(0, 0, 1);
	glNormal3dv(Normal(top[0], top[1], top[2]));
	for (int i = 0; i < 4; i++) {
		glVertex3dv(top[i].ToArray());
	}
	glEnd();
	glBegin(GL_QUADS);
	glColor3d(0, 1, 0);
	glNormal3dv(Normal(bot[0], bot[1], bot[2], true));
	for (int i = 0; i < 4; i++) {
		glVertex3dv(bot[i].ToArray());
	}
	glEnd();
}

Vec VecCompos(double* K1, double* K2) {
	vector<double> K = { 0, 0, 0 };
	double M[2][3];
	for (int i = 0; i < 3; i++) {
		M[0][i] = K1[i];
		M[1][i] = K2[i];
	}
	K[0] = (M[0][1] * M[1][2]) - (M[0][2] * M[1][1]);
	K[1] = ((M[0][0] * M[1][2]) - (M[0][2] * M[1][0])) * -1;
	K[2] = (M[0][0] * M[1][1]) - (M[0][1] * M[1][0]);
	return { K[0], K[1], K[2] };
}

double Cos(const Vec& K1, const Vec& K2) {
	return K1.Scal(K2) / (K1.Len() * K2.Len());
}

void _40(double delta, double t_max) {
	Vec P1 = { 0, 0, 0 };
	Vec P2 = { -4, 6, 7 };
	Vec P3 = { 4, 8, 2 };
	Vec P4 = { 10, 10, 0 };
	Beze(P1, P2, P3, P4);

	glPushMatrix();
	Vec point = F_Beze(P1, P2, P3, P4, t_max);
	Vec next = F_Beze(P1, P2, P3, P4, t_max + delta);

	Vec del = (next - point).Normaliz();

	Vec axisX = { 1, 0, 0 };
	Vec proecDelX = Vec(del.x, del.y, 0).Normaliz();
	double cosPsi = Cos(axisX, proecDelX);
	Vec vecSign = VecCompos(axisX.ToArray(), proecDelX.ToArray());
	double sign = vecSign.z / abs(vecSign.z);
	double psi = acos(cosPsi) * 180 / M_PI * sign;

	double teta = acos(del.z) * 180 / M_PI - 90;

	glTranslated(next.x, next.y, next.z);
	glRotated(psi, 0, 0, 1);
	glRotated(teta, 0, 1, 0);

	//Square({next.x, next.y, next.z});
	Square({ 0, 0, 0 });

	glPopMatrix();
}

array<Vec, 4> Triangle(double u, double v, double h) {
	array<Vec, 4> tmp;
	BezeSurfacePoint(u, v, tmp[0]);
	BezeSurfacePoint(u, v + h, tmp[1]);
	BezeSurfacePoint(u + h, v, tmp[2]);
	BezeSurfacePoint(u + h, v + h, tmp[3]);

	return tmp;
}

void Lines_Points() {
	glDisable(GL_LIGHTING);
	glColor3d(0.5, 0.3, 0.7);
	for (int i = 0; i < points.size(); i++) {
		glBegin(GL_LINE_STRIP);
		for (int j = 0; j < points[i].size(); j++)
			glVertex3dv(points[i][j].ToArray());
		glEnd();
	}	
	for (int i = 0; i < points.size(); i++) {
		glBegin(GL_LINE_STRIP);
		for (int j = 0; j < points[i].size(); j++)
			glVertex3dv(points[j][i].ToArray());
		glEnd();
	}
	glPointSize(15);
	glColor3d(1, 0, 0);
	glBegin(GL_POINTS);	
	for (const auto& i : points) {
		for(const auto& elem: i)
		glVertex3dv(elem.ToArray());
	}
	glEnd();
	glEnable(GL_LIGHTING);
}

void BezePlane() {
	double h = 0.05;
	for (double u = 0; u < 1.01- h; u += h) {
		glBegin(GL_TRIANGLES);
		for (double v = 0; v < 1.01 - h; v += h) {
			array<Vec, 4> arr = Triangle(u, v, h);

			glNormal3dv(Normal(arr[0], arr[1], arr[2]));
			glTexCoord3dv(TexPoint(arr[0]).ToArray());
			glVertex3dv(arr[0].ToArray());
			glTexCoord3dv(TexPoint(arr[1]).ToArray());
			glVertex3dv(arr[1].ToArray());
			glTexCoord3dv(TexPoint(arr[3]).ToArray());
			glVertex3dv(arr[3].ToArray());

			glNormal3dv(Normal(arr[0], arr[2], arr[3], true));
			glTexCoord3dv(TexPoint(arr[0]).ToArray());
			glVertex3dv(arr[0].ToArray());
			glTexCoord3dv(TexPoint(arr[2]).ToArray());
			glVertex3dv(arr[2].ToArray());
			glTexCoord3dv(TexPoint(arr[3]).ToArray());
			glVertex3dv(arr[3].ToArray());
		}
		glEnd();
	}
	Lines_Points();
}

double t_max = 0;
bool f = true;
double delta = 0.0018;
void Render(OpenGL* ogl)
{

	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);
	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);
	if (textureMode)
		glEnable(GL_TEXTURE_2D);

	if (lightMode)
		glEnable(GL_LIGHTING);


	//альфаналожение
	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


	//настройка материала
	GLfloat amb[] = { 0.2, 0.2, 0.1, 1. };
	GLfloat dif[] = { 0.4, 0.65, 0.5, 1. };
	GLfloat spec[] = { 0.9, 0.8, 0.3, 1. };
	GLfloat sh = 0.1f * 256;


	//фоновая
	glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
	//дифузная
	glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
	//зеркальная
	glMaterialfv(GL_FRONT, GL_SPECULAR, spec); \
		//размер блика
		glMaterialf(GL_FRONT, GL_SHININESS, sh);

	//чтоб было красиво, без квадратиков (сглаживание освещения)
	glShadeModel(GL_SMOOTH);
	//===================================
	//Прогать тут 
	glBindTexture(GL_TEXTURE_2D, texId);
	t_max += delta;
	if (t_max > 1) delta = -delta;
	if (t_max < 0)  delta = -delta;
	_40(delta, t_max);
	BezePlane();
	//===================================
	


	   //Сообщение вверху экрана


		glMatrixMode(GL_PROJECTION);	//Делаем активной матрицу проекций. 
										//(всек матричные операции, будут ее видоизменять.)
		glPushMatrix();   //сохраняем текущую матрицу проецирования (которая описывает перспективную проекцию) в стек 				    
		glLoadIdentity();	  //Загружаем единичную матрицу
		glOrtho(0, ogl->getWidth(), 0, ogl->getHeight(), 0, 1);	 //врубаем режим ортогональной проекции

		glMatrixMode(GL_MODELVIEW);		//переключаемся на модел-вью матрицу
		glPushMatrix();			  //сохраняем текущую матрицу в стек (положение камеры, фактически)
		glLoadIdentity();		  //сбрасываем ее в дефолт

		glDisable(GL_LIGHTING);



		GuiTextRectangle rec;		   //классик моего авторства для удобной работы с рендером текста.
		rec.setSize(300, 150);
		rec.setPosition(10, ogl->getHeight() - 150 - 10);


		std::stringstream ss;
		ss << "T - вкл/выкл текстур" << std::endl;
		ss << "L - вкл/выкл освещение" << std::endl;
		ss << "A - вкл/выкл альфаналожение" << std::endl;
		ss << "F - Свет из камеры" << std::endl;
		ss << "G - двигать свет по горизонтали" << std::endl;
		ss << "G+ЛКМ двигать свет по вертекали" << std::endl;
		ss << "Коорд. света: (" << light.pos.X() << ", " << light.pos.Y() << ", " << light.pos.Z() << ")" << std::endl;
		ss << "Коорд. камеры: (" << camera.pos.X() << ", " << camera.pos.Y() << ", " << camera.pos.Z() << ")" << std::endl;
		ss << "Параметры камеры: R=" << camera.camDist << ", fi1=" << camera.fi1 << ", fi2=" << camera.fi2 << std::endl;

		rec.setText(ss.str().c_str());
		rec.Draw();

		glMatrixMode(GL_PROJECTION);	  //восстанавливаем матрицы проекции и модел-вью обратьно из стека.
		glPopMatrix();


		glMatrixMode(GL_MODELVIEW);
		glPopMatrix();
}