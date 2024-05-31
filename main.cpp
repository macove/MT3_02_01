#include <Novice.h>
#include <cmath>
#include "Vector3.h"
#include <assert.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <imgui.h>

struct Matrix4x4
{
	float m[4][4];

};

struct Sphere {
	Vector3 center;
	float radius;

};

struct Segment {
	Vector3 origin;
	Vector3 diff;

};


Vector3 Add(const Vector3& v1, const Vector3& v2) {
	Vector3 result;

	result.x = v1.x + v2.x;
	result.y = v1.y + v2.y;
	result.z = v1.z + v2.z;

	return result;
}

Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result;

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			result.m[i][j] = 0;
			for (int k = 0; k < 4; ++k) {
				result.m[i][j] += m1.m[i][k] * m2.m[k][j];
			}
		}
	}
	return result;
}

Matrix4x4 MakeRotateXMatrix(float radian) {

	Matrix4x4 result;

	result = {
		1.0f,0.0f,0.0f,0.0f,
		0.0f,cosf(radian),sinf(radian),0.0f,
		0.0f,-sinf(radian),cosf(radian),0.0f,
		0.0f,0.0f,0.0f,1.0f };
	return result;
}

Matrix4x4 MakeRotateYMatrix(float radian) {

	Matrix4x4 result;

	result = {
		cosf(radian),0.0f,-sinf(radian),0.0f,
		0.0f,1.0f,0.0f,0.0f,
		sinf(radian),0.0f,cosf(radian),0.0f,
		0.0f,0.0f,0.0f,1.0f };
	return result;


}

Matrix4x4 MakeRotateZMatrix(float radian) {

	Matrix4x4 result;

	result = {
		cosf(radian),sinf(radian),0.0f,0.0f,
		-sinf(radian),cosf(radian),0.0f,0.0f,
		0.0f,0.0f,1.0f,0.0f,
		0.0f,0.0f,0.0f,1.0f };
	return result;


}

Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate) {

	Matrix4x4 rotateXMatrix = MakeRotateXMatrix(rotate.x);
	Matrix4x4 rotateYMatrix = MakeRotateYMatrix(rotate.y);
	Matrix4x4 rotateZMatrix = MakeRotateZMatrix(rotate.z);
	Matrix4x4 rotateXYZMatrix = Multiply(rotateXMatrix, Multiply(rotateYMatrix, rotateZMatrix));

	Matrix4x4 result;
	result = {
		scale.x * rotateXYZMatrix.m[0][0],scale.x * rotateXYZMatrix.m[0][1],scale.x * rotateXYZMatrix.m[0][2],0.0f,
		scale.y * rotateXYZMatrix.m[1][0],scale.y * rotateXYZMatrix.m[1][1],scale.y * rotateXYZMatrix.m[1][2],0.0f,
		scale.z * rotateXYZMatrix.m[2][0],scale.z * rotateXYZMatrix.m[2][1],scale.z * rotateXYZMatrix.m[2][2],0.0f,
		translate.x,translate.y,translate.z,1.0f
	};
	return result;
}

Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip) {
	Matrix4x4 result;

	result = {
		1 / aspectRatio * (1 / tanf(fovY / 2)),0.0f,0.0f,0.0f,
		0.0f,(1 / tanf(fovY / 2)),0.0f,0.0f,
		0.0f,0.0f,farClip / (farClip - nearClip),1.0f,
		0.0f,0.0f,-farClip * nearClip / (farClip - nearClip),0.0f
	};

	return result;

}

Matrix4x4 MakeOrthographicMatrix(float left, float top, float right, float bottom, float nearClip, float farClip) {

	Matrix4x4 result;
	result = {
	2 / (right - left),0.0f,0.0f,0.0f,
	0.0f,2 / (top - bottom),0.0f,0.0f,
	0.0f,0.0f,1 / (farClip - nearClip),0.0f,
	(left + right) / (left - right),(top + bottom) / (bottom - top),nearClip / (nearClip - farClip),1.0f
	};
	return result;
}

Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth) {

	Matrix4x4 result;

	result = {
	width / 2,0.0f,0.0f,0.0f,
	0.0f,-(height / 2),0.0f,0.0f,
	0.0f,0.0f,maxDepth - minDepth,0.0f,
	left + (width / 2),top + (height / 2),minDepth,1.0f
	};


	return result;
}

Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix) {

	float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] + vector.z * matrix.m[2][3] + 1.0f * matrix.m[3][3];
	assert(w != 0.0f);

	Vector3 result;
	result.x = vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] + vector.z * matrix.m[2][0] + 1.0f * matrix.m[3][0];
	result.y = vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] + vector.z * matrix.m[2][1] + 1.0f * matrix.m[3][1];
	result.z = vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] + vector.z * matrix.m[2][2] + 1.0f * matrix.m[3][2];

	result.x /= w;
	result.y /= w;
	result.z /= w;


	return result;

}

Matrix4x4 Inverse(const Matrix4x4& m) {
	Matrix4x4 result;
	float A;
	A = m.m[0][0] * m.m[1][1] * m.m[2][2] * m.m[3][3] + m.m[0][0] * m.m[1][2] * m.m[2][3] * m.m[3][1] + m.m[0][0] * m.m[1][3] * m.m[2][1] * m.m[3][2]
		- m.m[0][0] * m.m[1][3] * m.m[2][2] * m.m[3][1] - m.m[0][0] * m.m[1][2] * m.m[2][1] * m.m[3][3] - m.m[0][0] * m.m[1][1] * m.m[2][3] * m.m[3][2]
		- m.m[0][1] * m.m[1][0] * m.m[2][2] * m.m[3][3] - m.m[0][2] * m.m[1][0] * m.m[2][3] * m.m[3][1] - m.m[0][3] * m.m[1][0] * m.m[2][1] * m.m[3][2]
		+ m.m[0][3] * m.m[1][0] * m.m[2][2] * m.m[3][1] + m.m[0][2] * m.m[1][0] * m.m[2][1] * m.m[3][3] + m.m[0][1] * m.m[1][0] * m.m[2][3] * m.m[3][2]
		+ m.m[0][1] * m.m[1][2] * m.m[2][0] * m.m[3][3] + m.m[0][2] * m.m[1][3] * m.m[2][0] * m.m[3][1] + m.m[0][3] * m.m[1][1] * m.m[2][0] * m.m[3][2]
		- m.m[0][3] * m.m[1][2] * m.m[2][0] * m.m[3][1] - m.m[0][2] * m.m[1][1] * m.m[2][0] * m.m[3][3] - m.m[0][1] * m.m[1][3] * m.m[2][0] * m.m[3][2]
		- m.m[0][1] * m.m[1][2] * m.m[2][3] * m.m[3][0] - m.m[0][2] * m.m[1][3] * m.m[2][1] * m.m[3][0] - m.m[0][3] * m.m[1][1] * m.m[2][2] * m.m[3][0]
		+ m.m[0][3] * m.m[1][2] * m.m[2][1] * m.m[3][0] + m.m[0][2] * m.m[1][1] * m.m[2][3] * m.m[3][0] + m.m[0][1] * m.m[1][3] * m.m[2][2] * m.m[3][0];

	result.m[0][0] = 1 / A * (m.m[1][1] * m.m[2][2] * m.m[3][3] + m.m[1][2] * m.m[2][3] * m.m[3][1] + m.m[1][3] * m.m[2][1] * m.m[3][2]
		- m.m[1][3] * m.m[2][2] * m.m[3][1] - m.m[1][2] * m.m[2][1] * m.m[3][3] - m.m[1][1] * m.m[2][3] * m.m[3][2]);
	result.m[0][1] = 1 / A * (-m.m[0][1] * m.m[2][2] * m.m[3][3] - m.m[0][2] * m.m[2][3] * m.m[3][1] - m.m[0][3] * m.m[2][1] * m.m[3][2]
		+ m.m[0][3] * m.m[2][2] * m.m[3][1] + m.m[0][2] * m.m[2][1] * m.m[3][3] + m.m[0][1] * m.m[2][3] * m.m[3][2]);
	result.m[0][2] = 1 / A * (m.m[0][1] * m.m[1][2] * m.m[3][3] + m.m[0][2] * m.m[1][3] * m.m[3][1] + m.m[0][3] * m.m[1][1] * m.m[3][2]
		- m.m[0][3] * m.m[1][2] * m.m[3][1] - m.m[0][2] * m.m[1][1] * m.m[3][3] - m.m[0][1] * m.m[1][3] * m.m[3][2]);
	result.m[0][3] = 1 / A * (-m.m[0][1] * m.m[1][2] * m.m[2][3] - m.m[0][2] * m.m[1][3] * m.m[2][1] - m.m[0][3] * m.m[1][1] * m.m[2][2]
		+ m.m[0][3] * m.m[1][2] * m.m[2][1] + m.m[0][2] * m.m[1][1] * m.m[2][3] + m.m[0][1] * m.m[1][3] * m.m[2][2]);

	result.m[1][0] = 1 / A * (-m.m[1][0] * m.m[2][2] * m.m[3][3] - m.m[1][2] * m.m[2][3] * m.m[3][0] - m.m[1][3] * m.m[2][0] * m.m[3][2]
		+ m.m[1][3] * m.m[2][2] * m.m[3][0] + m.m[1][2] * m.m[2][0] * m.m[3][3] + m.m[1][0] * m.m[2][3] * m.m[3][2]);
	result.m[1][1] = 1 / A * (m.m[0][0] * m.m[2][2] * m.m[3][3] + m.m[0][2] * m.m[2][3] * m.m[3][0] + m.m[0][3] * m.m[2][0] * m.m[3][2]
		- m.m[0][3] * m.m[2][2] * m.m[3][0] - m.m[0][2] * m.m[2][0] * m.m[3][3] - m.m[0][0] * m.m[2][3] * m.m[3][2]);
	result.m[1][2] = 1 / A * (-m.m[0][0] * m.m[1][2] * m.m[3][3] - m.m[0][2] * m.m[1][3] * m.m[3][0] - m.m[0][3] * m.m[1][0] * m.m[3][2]
		+ m.m[0][3] * m.m[1][2] * m.m[3][0] + m.m[0][2] * m.m[1][0] * m.m[3][3] + m.m[0][0] * m.m[1][3] * m.m[3][2]);
	result.m[1][3] = 1 / A * (m.m[0][0] * m.m[1][2] * m.m[2][3] + m.m[0][2] * m.m[1][3] * m.m[2][0] + m.m[0][3] * m.m[1][0] * m.m[2][2]
		- m.m[0][3] * m.m[1][2] * m.m[2][0] - m.m[0][2] * m.m[1][0] * m.m[2][3] - m.m[0][0] * m.m[1][3] * m.m[2][2]);

	result.m[2][0] = 1 / A * (m.m[1][0] * m.m[2][1] * m.m[3][3] + m.m[1][1] * m.m[2][3] * m.m[3][0] + m.m[1][3] * m.m[2][0] * m.m[3][1]
		- m.m[1][3] * m.m[2][1] * m.m[3][0] - m.m[1][1] * m.m[2][0] * m.m[3][3] - m.m[1][0] * m.m[2][3] * m.m[3][1]);
	result.m[2][1] = 1 / A * (-m.m[0][0] * m.m[2][1] * m.m[3][3] - m.m[0][1] * m.m[2][3] * m.m[3][0] - m.m[0][3] * m.m[2][0] * m.m[3][1]
		+ m.m[0][3] * m.m[2][1] * m.m[3][0] + m.m[0][1] * m.m[2][0] * m.m[3][3] + m.m[0][0] * m.m[2][3] * m.m[3][1]);
	result.m[2][2] = 1 / A * (m.m[0][0] * m.m[1][1] * m.m[3][3] + m.m[0][1] * m.m[1][3] * m.m[3][0] + m.m[0][3] * m.m[1][0] * m.m[3][1]
		- m.m[0][3] * m.m[1][1] * m.m[3][0] - m.m[0][1] * m.m[1][0] * m.m[3][3] - m.m[0][0] * m.m[1][3] * m.m[3][1]);
	result.m[2][3] = 1 / A * (-m.m[0][0] * m.m[1][1] * m.m[2][3] - m.m[0][1] * m.m[1][3] * m.m[2][0] - m.m[0][3] * m.m[1][0] * m.m[2][1]
		+ m.m[0][3] * m.m[1][1] * m.m[2][0] + m.m[0][1] * m.m[1][0] * m.m[2][3] + m.m[0][0] * m.m[1][3] * m.m[2][1]);

	result.m[3][0] = 1 / A * (-m.m[1][0] * m.m[2][1] * m.m[3][2] - m.m[1][1] * m.m[2][2] * m.m[3][0] - m.m[1][2] * m.m[2][0] * m.m[3][1]
		+ m.m[1][2] * m.m[2][1] * m.m[3][0] + m.m[1][1] * m.m[2][0] * m.m[3][2] + m.m[1][0] * m.m[2][2] * m.m[3][1]);
	result.m[3][1] = 1 / A * (m.m[0][0] * m.m[2][1] * m.m[3][2] + m.m[0][1] * m.m[2][2] * m.m[3][0] + m.m[0][2] * m.m[2][0] * m.m[3][1]
		- m.m[0][2] * m.m[2][1] * m.m[3][0] - m.m[0][1] * m.m[2][0] * m.m[3][2] - m.m[0][0] * m.m[2][2] * m.m[3][1]);
	result.m[3][2] = 1 / A * (-m.m[0][0] * m.m[1][1] * m.m[3][2] - m.m[0][1] * m.m[1][2] * m.m[3][0] - m.m[0][2] * m.m[1][0] * m.m[3][1]
		+ m.m[0][2] * m.m[1][1] * m.m[3][0] + m.m[0][1] * m.m[1][0] * m.m[3][2] + m.m[0][0] * m.m[1][2] * m.m[3][1]);
	result.m[3][3] = 1 / A * (m.m[0][0] * m.m[1][1] * m.m[2][2] + m.m[0][1] * m.m[1][2] * m.m[2][0] + m.m[0][2] * m.m[1][0] * m.m[2][1]
		- m.m[0][2] * m.m[1][1] * m.m[2][0] - m.m[0][1] * m.m[1][0] * m.m[2][2] - m.m[0][0] * m.m[1][2] * m.m[2][1]);
	return result;

}

void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix) {
	const float kGridHalfWidth = 2.0f; //Gridの半分の幅
	const uint32_t kSubdivision = 10;  //分割数
	const float kGridEvery = (kGridHalfWidth * 2.0f) / float(kSubdivision); //1つ分の長さ
	for (uint32_t xIndex = 0; xIndex <= kSubdivision; ++xIndex) {
		//上の情報を使ってワールド座標系上の始点と終点を求める
		//スクリーン座標系まで変換をかける
		//変換した座標を使って表示	

		Vector3 startPoint(kGridHalfWidth - (kGridEvery * xIndex), 0.0f, -kGridHalfWidth);
		Vector3 endPoint(kGridHalfWidth - (kGridEvery * xIndex), 0, kGridHalfWidth);

		Vector3 screenStartPoint = Transform(startPoint, viewProjectionMatrix);
		Vector3 screenEndPoint = Transform(endPoint, viewProjectionMatrix);

		Vector3 screenStartPoint1 = Transform(screenStartPoint, viewportMatrix);
		Vector3 screenEndPoint1 = Transform(screenEndPoint, viewportMatrix);

		Novice::DrawLine(int(screenStartPoint1.x), int(screenStartPoint1.y), int(screenEndPoint1.x), int(screenEndPoint1.y), 0xAAAAAAFF);


	}
	for (uint32_t zIndex = 0; zIndex <= kSubdivision; ++zIndex) {

		Vector3 startPoint(-kGridHalfWidth, 0.0f, kGridHalfWidth - (kGridEvery * zIndex));
		Vector3 endPoint(kGridHalfWidth, 0, kGridHalfWidth - (kGridEvery * zIndex));

		Vector3 screenStartPoint = Transform(startPoint, viewProjectionMatrix);
		Vector3 screenEndPoint = Transform(endPoint, viewProjectionMatrix);

		Vector3 screenStartPoint1 = Transform(screenStartPoint, viewportMatrix);
		Vector3 screenEndPoint1 = Transform(screenEndPoint, viewportMatrix);

		Novice::DrawLine(int(screenStartPoint1.x), int(screenStartPoint1.y), int(screenEndPoint1.x), int(screenEndPoint1.y), 0xAAAAAAFF);

	}
}
Vector3 Subtract(const Vector3& v1, const Vector3& v2) {
	Vector3 result;

	result.x = v1.x - v2.x;
	result.y = v1.y - v2.y;
	result.z = v1.z - v2.z;

	return result;
}
Vector3 Project(const Vector3& v1, const Vector3& v2)
{
	Vector3 result;
	float a = ((v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z));
	float b = ((v2.x * v2.x) + (v2.y * v2.y) + (v2.z * v2.z));
	result.x = a / b * v2.x;
	result.y = a / b * v2.y;
	result.z = a / b * v2.z;

	return result;
}
Vector3 ClosestPoint(const Vector3& point, const Segment& segment) {

	Vector3 project = Project(Subtract(point, segment.origin), segment.diff);

	project = Add(project, segment.origin);

	return project;
}

void DrawSphere(const Sphere& sphere, const Matrix4x4& viewProjection, const Matrix4x4 viewportMatrix, uint32_t color) {
	const uint32_t kSubdivision = 16;
	const float kLonEvery = float(M_PI * 2) / kSubdivision;
	const float kLatEvery = float(M_PI) / kSubdivision;
	for (uint32_t latIndex = 0; latIndex < kSubdivision; ++latIndex)
	{
		float lat = float(-M_PI / 2.0f) + (kLatEvery * latIndex);

		for (uint32_t lonIndex = 0; lonIndex < kSubdivision; ++lonIndex)
		{
			float lon = float(lonIndex * kLonEvery);

			Vector3 a, b, c;

			a = { (cosf(lat) * cosf(lon) * sphere.radius + sphere.center.x),(sinf(lat) * sphere.radius + sphere.center.y) ,(cosf(lat) * sinf(lon) * sphere.radius + sphere.center.z) };
			b = { (cosf(lat + kLatEvery) * cosf(lon) * sphere.radius + sphere.center.x) ,(sinf(lat + kLatEvery) * sphere.radius + sphere.center.y) ,(cosf(lat + kLatEvery) * sinf(lon) * sphere.radius + sphere.center.z) };
			c = { (cosf(lat) * cosf(lon + kLonEvery) * sphere.radius + sphere.center.x) ,(sinf(lat) * sphere.radius + sphere.center.y) ,(cosf(lat) * sinf(lon + kLonEvery) * sphere.radius + sphere.center.z) };

			Vector3 VecA = Transform(Transform(a, viewProjection), viewportMatrix);
			Vector3 VecB = Transform(Transform(b, viewProjection), viewportMatrix);
			Vector3 VecC = Transform(Transform(c, viewProjection), viewportMatrix);

			Novice::DrawLine(int(VecA.x), int(VecA.y), int(VecB.x), int(VecB.y), color);
			Novice::DrawLine(int(VecA.x), int(VecA.y), int(VecC.x), int(VecC.y), color);

		}
	}

}

bool IsCollision(const Sphere& s1, const Sphere& s2) {

	float a = s1.center.x - s2.center.x;
	float b = s1.center.y - s2.center.y;
	float c = s1.center.z - s2.center.z;
	float distance = sqrtf(a * a + b * b + c * c);

	if (distance <= s1.radius + s2.radius) {
		return true;
	}
	else {
		return false;
	}

}

const char kWindowTitle[] = "学籍番号";

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	// キー入力結果を受け取る箱
	char keys[256] = { 0 };
	char preKeys[256] = { 0 };

	Vector3 rotate{ 0.0f,0.0f,0.0f };

	Vector3 translate{ 0.0f,0.0f,0.0f };

	Vector3 cameraPosition{ 0.0f,1.9f,-6.49f };

	Vector3 cameraRotate{ 0.26f,0.0f,0.0f };

	Segment segment{ {-2.0f,-1.0f,0.0f},{3.0f,2.0f,2.0f} };
	Vector3 point{ -1.5f,0.6f,0.6f };

	Sphere A{ 0.0f,0.0f,0.0f,0.5f };
	Sphere B{ 1.0f,1.0f,0.0f,0.3f };

	int color = 0xFFFFFFFF;

	int kWindowWidth = 1280;
	int kWindowHeight = 720;

	//float cameraRadius = 10.0f; 
	//float cameraTheta = float(M_PI / 36); 
	//float cameraPhi = float(M_PI / 36); 

	
	
	

	//Vector3 cameraPositionA{ 0.0f,1.9f,-6.49f };

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///
		
		//cameraPositionA.x = cameraRadius * sinf(cameraTheta) * cosf(cameraPhi);
		//cameraPositionA.y = cameraRadius * sinf(cameraTheta) * sinf(cameraPhi);
		//cameraPositionA.z = cameraRadius * cosf(cameraTheta);

		Matrix4x4 worldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, rotate, translate);
		Matrix4x4 cameraMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, cameraRotate, cameraPosition);
		Matrix4x4 viewMatrix = Inverse(cameraMatrix);
		Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
		Matrix4x4 worldViewProjectionMatrix = Multiply(worldMatrix, Multiply(viewMatrix, projectionMatrix));
		Matrix4x4 viewProjectionMatrix = Multiply(viewMatrix, projectionMatrix);
		Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);


		if (IsCollision(A, B) == true) {

			color = 0xFF0000FF;
		} else {
			color = 0xFFFFFFFF;
		}





		//Matrix4x4 cameraMatrixA = MakeAffineMatrix({ 1.0f, 1.0f, 1.0f }, cameraRotate, cameraPositionA);
		//Matrix4x4 viewMatrixA = Inverse(cameraMatrixA);

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///
		DrawGrid(viewProjectionMatrix, viewportMatrix);
		DrawSphere(A, viewProjectionMatrix, viewportMatrix, color);
		DrawSphere(B, viewProjectionMatrix, viewportMatrix, 0xFFFFFFFF);
		


		ImGui::Begin("Window");
		//ImGui::DragFloat3("cameraPositionA", &cameraPositionA.x, 0.03f);
		ImGui::DragFloat3("positionA", &A.center.x, 0.03f);
		ImGui::DragFloat("A radius", &A.radius, 0.01f);
		ImGui::DragFloat3("positionB", &B.center.x, 0.03f);
		ImGui::DragFloat("B radius", &B.radius, 0.01f);
		ImGui::End();

		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}
