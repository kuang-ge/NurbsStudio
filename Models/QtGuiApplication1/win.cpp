#include "win.h"
#include <QDebug>
#include <QTimer>
#include <QKeyEvent>
#include <QDateTime>
#include "MyDoc.h"
#include "Fragment.h"
#include "qapplication.h"

win::win(QWidget *parent) : QOpenGLWidget(parent)
{
	MyDoc::Ptr pdoc = MyDoc::getInstance();
	//setGeometry(0, 0, 800, 800);
	camera = std::make_unique<Camera>(QVector3D(pdoc->cameraX, pdoc->cameraY, pdoc->cameraZ));//������ĳ�ʼλ��
	m_bLeftPressed = false;
	QSurfaceFormat surfaceFormat;
	surfaceFormat.setSamples(6);//���ز���
	setFormat(surfaceFormat); //setFormat��QOpenGLWidget�ĺ���
	m_pTimer = new QTimer(this);
	connect(m_pTimer, &QTimer::timeout, this, [=] {
		m_nTimeValue += 1;
		update();
	});
	m_pTimer->start(10);//����������fps����ֵԽ��fpsԽ��
}

win::~win()
{
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);
}

void win::initializeGL() {

	this->initializeOpenGLFunctions();
	createShader();

	//�Ȱ�����ϵ���ú�
	//----------
	float X[] = {

		0.0f,0.0f,0.0f,  1.0f,0.0f,0.0f,
		1000.0f,0.0f,0.0f,  1.0f,0.0f,0.0f,

		0.0f,0.0f,0.0f,  0.0f,1.0f,0.0f,
		0.0f,1000.0f,0.0f,  0.0f,1.0f,0.0f,

		0.0f,0.0f,0.0f,  0.0f,0.0f,1.0f,
		0.0f,0.0f,1000.0f,  0.0f,0.0f,1.0f

	};
	//setGeometry(300, 300, 800, 800);

	//������Ķ����������
	//---------
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(X), X, GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)0);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	
	// ����OpenGL��ȫ��״̬
	// --------
	glEnable(GL_DEPTH_TEST);//������Ȼ���
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void win::resizeGL(int w, int h) {
	glViewport(0,0 ,w, h);
}

void win::paintGL() {

	//�洢�ļ��Ļ�ȡ
	// -----
	MyDoc::Ptr pdoc = MyDoc::getInstance();
	// �趨ÿ�μ��̰���֮����ٶ�
	// -----
	camera->processInput(pdoc->_keyBordSensitive);//�ٶ�
	camera->setMouseSpeed(pdoc->_mouseSenstive);
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//�������������λ��
	//----------
	if (pdoc->m_restView == true)
	{
		camera->setLocation(QVector3D(pdoc->cameraX, pdoc->cameraY, pdoc->cameraZ));
		pdoc->m_restView = false;
	}

	//��ʾģʽѡ��
	//----------
	switch (pdoc->m_GLMode)//ѡ����ʾ��ģʽ���ߣ��棬��
	{
	case 0:
		glPointSize(pdoc->m_dotSZ);
		glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
		break;
	case 3:
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		break;
	case 7:
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		break;
	default:
		break;
	}

	//������������
	//----------
	triangleShader.bind();
	QMatrix4x4 projection;
	projection.perspective(camera->zoom, 1.0f * width() / height(), 0.1f, 2000.0f);
	triangleShader.setUniformValue("projection", projection);
	QMatrix4x4 view;
	view = camera->getViewMatrix();
	triangleShader.setUniformValue("view", view);
	QMatrix4x4 model;
	//����ctrl��֮������������ת
	model.rotate(-ro.x(), 1, 0, 0);
	model.rotate(ro.y(), 0, 1, 0);
	model.rotate(ro.z(), 0, 0, 1);
	//��ת����
	triangleShader.setUniformValue("model", model);
	triangleShader.setUniformValue("ourColor", pdoc->_surfaceColor);

	//�Ƿ���Ҫȥ��Ƭ��ʾ
	//----------
	vector<QVector3D> n = { {0.3,0.5,0.5},{0.5,1.0,0.5},{0.5,0.5,1.0},{0.5,1.0,1.0},
	{1.0,0.5,1.0}, {0.5,0.8,0.5} };
	if (pdoc->m_CellIdx <= -1 && pdoc->m_Fragment == false)//��Ҫ��Ƭ��ʾ
	{
		for (int i = 0; i != pdoc->VAO_D_ARRAY.size(); ++i)
		{
			triangleShader.setUniformValue("ourColor", n[i%6]);
			glBindVertexArray(pdoc->VAO_D_ARRAY[i]);
			glDrawElements(GL_TRIANGLES, pdoc->EACH_VAO_D_POINT_NUMBER[i], GL_UNSIGNED_INT, 0);
		}
	}
	else if (pdoc->m_CellIdx >= 0)//Ҫ��Ƭ��ʾ
	{
		for (int i = pdoc->m_CellIdx * 6; i != pdoc->m_CellIdx * 6 + 6; ++i)
		{
			glBindVertexArray(pdoc->VAO_D_ARRAY[i]);
			glDrawElements(GL_TRIANGLES, pdoc->EACH_VAO_D_POINT_NUMBER[i], GL_UNSIGNED_INT, 0);
		}
	}

	//��������
	//----------
	triangleShader.setUniformValue("ourColor", pdoc->_bordLineColor);
	glLineWidth(pdoc->_lineWidth);
	for (int i = 0; i != pdoc->VAO_L_ARRAY.size(); ++i)
	{
		glBindVertexArray(pdoc->VAO_L_ARRAY[i]);
		glDrawArrays(GL_LINE_STRIP, 0, pdoc->EACH_VAO_L_POINT_NUMBER[i]);
	}

	//�����Ƶ�
	//----------
	if (pdoc->m_showCtrlPts == true)//���������ʾ���Ƶ㣬���Ǿ���ʾ���Ƶ�
	{
		triangleShader.setUniformValue("ourColor", pdoc->_controlPoinrColor);
		glPointSize(pdoc->m_CPTSdotSZ);
		if (pdoc->m_CellIdx == -1)//��ʾȫ���Ŀ��Ƶ�
		{
			for (int i = 0; i != pdoc->VAO_CP_ARRAY.size(); ++i)
			{
				glBindVertexArray(pdoc->VAO_CP_ARRAY[i]);
				glDrawArrays(GL_POINTS, 0, pdoc->EACH_VAO_CP_POINT_NUMBER[i]);
			}
		}
		else//��ʾѡ����ģ�͵Ŀ��Ƶ�
		{
			for (int i = pdoc->m_CellIdx /** 6*/; i != pdoc->m_CellIdx+1 /** 6 + 6*/; ++i)
			{
				glBindVertexArray(pdoc->VAO_CP_ARRAY[i]);
				glDrawArrays(GL_POINTS, 0, pdoc->EACH_VAO_CP_POINT_NUMBER[i]);
			}
		}

		
	}

	//�����Ҫ��Ƭ
	//----------
	if (pdoc->m_Fragment == true)
	{
		glPointSize(10);
		triangleShader.setUniformValue("ourColor", pdoc->_sliceSurfaceColor);
		
		for (int i = 0; i != pdoc->VAO_FRAGMENT_ARRAY.size(); ++i)
		{
			glBindVertexArray(pdoc->VAO_FRAGMENT_ARRAY[i]);
			glDrawElements(GL_TRIANGLES, pdoc->EACH_VAO_F_NUMBER[i], GL_UNSIGNED_INT, 0);
		}

		triangleShader.setUniformValue("ourColor", pdoc->_sliceSuppotrColor);
		for (int i = 0; i != pdoc->VAO_SUPPORT_ARRAY.size(); ++i)
		{
			glBindVertexArray(pdoc->VAO_SUPPORT_ARRAY[i]);
			glDrawArrays(GL_LINES, 0, pdoc->EACH_VAO_SUPPORT_NUMBER[i]);
		}
	}
	triangleShader.release();

	//������ϵ
	//----------
	coorShader.bind();
	if (pdoc->m_showCoordinates == true)
	{
		coorShader.setUniformValue("projection", projection);
		coorShader.setUniformValue("view", view);
		QMatrix4x4 model1;
		coorShader.setUniformValue("model", model1);
		glLineWidth(1.0f);
		glBindVertexArray(VAO);
		glDrawArrays(GL_LINE_STRIP, 0, 6);
		glDrawArrays(GL_POINTS, 0, 100000);

		//�������������ɫ�����зǾ���ģ�͵���Ⱦ
		/*glPointSize(10.0f);
		for (int i = 0; i != pdoc->VAO_D.size(); ++i)
		{
			glBindVertexArray(pdoc->VAO_D[i]);
			glDrawElements(GL_TRIANGLES, pdoc->VAO_D_SIZE[i], GL_UNSIGNED_INT, 0);
		}*/
	}
	coorShader.release();
}

void win::keyPressEvent(QKeyEvent *event)
{
	int key = event->key();
	if (key >= 0 && key < 1024)
		camera->keys[key] = true;
}

void win::keyReleaseEvent(QKeyEvent *event)
{
	int key = event->key();
	if (key >= 0 && key < 1024)
		camera->keys[key] = false;
}

void win::mousePressEvent(QMouseEvent *event)
{
	if (event->button() == Qt::LeftButton) {
		m_bLeftPressed = true;
		m_lastPos = event->pos();
	}
}

void win::mouseReleaseEvent(QMouseEvent *event)
{
	Q_UNUSED(event);

	m_bLeftPressed = false;
}

void win::mouseMoveEvent(QMouseEvent *event)
{
	int xpos = event->pos().x();
	int ypos = event->pos().y();

	int xoffset = xpos - m_lastPos.x();
	int yoffset = m_lastPos.y() - ypos;
	m_lastPos = event->pos();

	if ((QApplication::keyboardModifiers() == Qt::ControlModifier))
	{
		float a = ro.x() + yoffset / 10.0;
		ro.setX(a);
		if (ro.x() > 360.0 || ro.x() < -360.0)
			ro.setX(0.0);
		float b = ro.y() + xoffset / 10.0;
		ro.setY(b);
		if (ro.y() > 360.0 || ro.y() < -360.0)
			ro.setY(0.0);
	}
	else
	{
		camera->processMouseMovement(xoffset, yoffset);
	}
}

void win::wheelEvent(QWheelEvent *event)
{
	QPoint offset = event->angleDelta();
	camera->processMouseScroll(offset.y() / 20.0f);
}

bool win::createShader()
{
	
	bool success = triangleShader.addShaderFromSourceFile(QOpenGLShader::Vertex, ":/triangle.vert");
	if (!success) {
		qDebug() << "shaderProgram addShaderFromSourceFile failed!" << triangleShader.log();
		return success;
	}

	success = triangleShader.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/triangle.frag");
	if (!success) {
		qDebug() << "shaderProgram addShaderFromSourceFile failed!" << triangleShader.log();
		return success;
	}

	success = triangleShader.link();
	if (!success) {
		qDebug() << "shaderProgram link failed!" << triangleShader.log();
	}

	success = coorShader.addShaderFromSourceFile(QOpenGLShader::Vertex, ":/coor.vert");
	if (!success) {
		qDebug() << "shaderProgram addShaderFromSourceFile failed!" << coorShader.log();
		return success;
	}

	success = coorShader.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/coor.frag");
	if (!success) {
		qDebug() << "shaderProgram addShaderFromSourceFile failed!" << coorShader.log();
		return success;
	}

	success = coorShader.link();
	if (!success) {
		qDebug() << "shaderProgram link failed!" << coorShader.log();
	}

	return success;
}

uint win::loadTexture(const QString& path)
{
	uint textureID;
	glGenTextures(1, &textureID);

	QImage image = QImage(path).convertToFormat(QImage::Format_RGBA8888).mirrored(true, true);
	if (!image.isNull()) {
		glBindTexture(GL_TEXTURE_2D, textureID);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, image.width(), image.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, image.bits());
		glGenerateMipmap(GL_TEXTURE_2D);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	}

	return textureID;
}

//void win::createData(QOpenGLContext * p)
//{
//	//���õ�ǰ����Ⱦ����
//	QSurface *a = p->surface();
//	p->makeCurrent(a);
//
//	MyDoc::Ptr pdoc = MyDoc::getInstance();
//	pdoc->m_CellIdx = -1;
//
//	std::vector<unsigned int> indices;
//	std::vector<unsigned int> base = { 0,1,2,2,3,0 };
//
//	//�������Ŵ���һ��һά���飬����ԭ��������
//	//��������ݱ���
//	for (int i = 0; i != pdoc->m_ShowData.size(); ++i)
//	{
//		int pointsize = 0;
//		for (int j = 0; j != pdoc->m_ShowData[i].size(); ++j)
//		{
//			for (int k = 0; k != pdoc->m_ShowData[i][j].size(); ++k)
//			{
//				pdoc->showDta_vertices.push_back(pdoc->m_ShowData[i][j][k].x);
//				pdoc->showDta_vertices.push_back(pdoc->m_ShowData[i][j][k].y);
//				pdoc->showDta_vertices.push_back(pdoc->m_ShowData[i][j][k].z);
//			}
//			for (int k = 0; k != base.size(); ++k)
//			{
//				indices.push_back(base[k]);
//				base[k] += 4;
//			}
//			pointsize += 6;
//		}
//
//		uint VAO, VBO, EBO;
//		glGenVertexArrays(1, &VAO);//�����������
//		glGenBuffers(1, &VBO);//��������
//		glGenBuffers(1, &EBO);//��������
//		glBindVertexArray(VAO);//��ʼ��
//		glBindBuffer(GL_ARRAY_BUFFER, VBO);
//		glBufferData(GL_ARRAY_BUFFER, pdoc->showDta_vertices.size() * sizeof(float), pdoc->showDta_vertices.data(), GL_STATIC_DRAW);
//		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
//		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(float)*indices.size(), indices.data(), GL_STATIC_DRAW);
//		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid *)0);
//		glEnableVertexAttribArray(0);
//		glBindBuffer(GL_ARRAY_BUFFER, 0);
//		glBindVertexArray(0);
//		pdoc->EACH_VAO_D_POINT_NUMBER.push_back(pointsize);
//		pdoc->VAO_D_ARRAY.push_back(VAO);
//		pdoc->showDta_vertices.clear();
//		indices.clear();
//		base = { 0,1,2,2,3,0 };
//	}
//
//	//�����ߵ����ݣ��߲���Ҫ�������飬����Ϊ������Ҫһ��һ��ķֱ𻭳���
//	for (int i = 0; i != pdoc->m_ShowLines.size(); ++i)
//	{
//		for (int j = 0; j != pdoc->m_ShowLines[i].size(); ++j)
//		{
//			int pointSize = 0;//������һ���߶������ж��ٸ���
//			for (int k = 0; k != pdoc->m_ShowLines[i][j].size(); ++k)
//			{
//				pdoc->showLines_vertices.push_back(pdoc->m_ShowLines[i][j][k].x);
//				pdoc->showLines_vertices.push_back(pdoc->m_ShowLines[i][j][k].y);
//				pdoc->showLines_vertices.push_back(pdoc->m_ShowLines[i][j][k].z);
//				pointSize++;
//			}
//			unsigned int VBO, VAO;
//			glGenVertexArrays(1, &VAO);
//			glGenBuffers(1, &VBO);
//			glBindVertexArray(VAO);
//			glBindBuffer(GL_ARRAY_BUFFER, VBO);
//			glBufferData(GL_ARRAY_BUFFER, pdoc->showLines_vertices.size() * sizeof(float), pdoc->showLines_vertices.data(), GL_STATIC_DRAW);
//			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid *)0);
//			glEnableVertexAttribArray(0);
//			glBindBuffer(GL_ARRAY_BUFFER, 0);
//			glBindVertexArray(0);
//			pdoc->VAO_L_ARRAY.push_back(VAO);//����һ��Ķ��㻺������¼����
//			pdoc->EACH_VAO_L_POINT_NUMBER.push_back(pointSize);//����һ���ж��ٸ����¼����
//			pdoc->showLines_vertices.clear();//������գ�����������һ��ѭ��
//		}
//	}
//	//���Ƶ������
//	for (int i = 0; i != pdoc->m_ShowData4D.size(); ++i)
//	{
//		int pointsize = 0;
//		for (int j = 0; j != pdoc->m_ShowData4D[i].size(); ++j)
//		{
//			pdoc->showCp_vertices.push_back(pdoc->m_ShowData4D[i][j].x);
//			pdoc->showCp_vertices.push_back(pdoc->m_ShowData4D[i][j].y);
//			pdoc->showCp_vertices.push_back(pdoc->m_ShowData4D[i][j].z);
//			pointsize++;
//		}
//		unsigned int VAO, VBO;
//		glGenVertexArrays(1, &VAO);
//		glGenBuffers(1, &VBO);
//		glBindVertexArray(VAO);
//		glBindBuffer(GL_ARRAY_BUFFER, VBO);
//		glBufferData(GL_ARRAY_BUFFER, pdoc->showCp_vertices.size() * sizeof(float), pdoc->showCp_vertices.data(), GL_STATIC_DRAW);
//		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid *)0);
//		glEnableVertexAttribArray(0);
//		glBindBuffer(GL_ARRAY_BUFFER, 0);
//		glBindVertexArray(0);
//		pdoc->EACH_VAO_CP_POINT_NUMBER.push_back(pointsize);
//		pdoc->VAO_CP_ARRAY.push_back(VAO);
//		pdoc->showCp_vertices.clear();
//	}
//}

void win::createFragment()
{
	makeCurrent();//�趨��ǰ����Ⱦ����

	MyDoc::Ptr pdoc = MyDoc::getInstance();
	pdoc->m_Fragment = !pdoc->m_Fragment;

	//����������һ����Ƭ�࣬ƽ����xoy�棬��begin��ʼ��end�������м��з�degree��ƽ��
	Fragment f(pdoc->m_fragmentA, pdoc->m_fragmentB, pdoc->m_fragmentC, pdoc->m_fragmentBegin, pdoc->m_fragmentEnd, pdoc->m_fragmentDegree);

	vector<vector<float>> temp;
	temp = f.getPoints(pdoc->m_ShowData);//���ܷ��صĶ������飬�ֱ���ÿ����������е�
	for (int i = 0; i != temp.size(); ++i)//���ǰ�ÿ������ֿ���ʾ
	{
		vector<int> indices;//�������һ���������飬���㻭ͼ
		vector<int> Indices = { 0,1,2 };
		for (int j = 0; j != temp[i].size() / 3 - 2; ++j)
		{
			indices.insert(indices.end(), Indices.begin(), Indices.end());
			Indices[1] += 1;
			Indices[2] += 1;
		}
		unsigned int VAO, VBO, EBO;
		glGenVertexArrays(1, &VAO);//�����������
		glGenBuffers(1, &VBO);//��������
		glGenBuffers(1, &EBO);//��������
		glBindVertexArray(VAO);//��ʼ��
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, temp[i].size() * sizeof(float), temp[i].data(), GL_STATIC_DRAW);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(float)*indices.size(), indices.data(), GL_STATIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid *)0);
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);
		pdoc->VAO_FRAGMENT_ARRAY.push_back(VAO);
		pdoc->EACH_VAO_F_NUMBER.push_back(temp[i].size());
	}

	vector<float> aa = f.getSupportPoints();
	unsigned int VAO, VBO, EBO;
	glGenVertexArrays(1, &VAO);//�����������
	glGenBuffers(1, &VBO);//��������
	glBindVertexArray(VAO);//��ʼ��
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, aa.size() * sizeof(float), aa.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid *)0);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
	pdoc->VAO_SUPPORT_ARRAY.push_back(VAO);
	pdoc->EACH_VAO_SUPPORT_NUMBER.push_back(aa.size() / 3);
}