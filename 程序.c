//��������еԭ��4-12�����������4-1��ֱ���Ӷ�������͹����Ƴ���
//֧�֣�utf-8��c��windows
//���ߣ���ϭ��
//build��2020/05/12
//updata��2020/05/13
//************************************
//�汾��v1.0
//��������
//�����������
//�󵼺�������bug
//************************************
//v1.1
//�޸����󵼺�����bug
//�޸��˷ֶκ�����bug
//ʵ���������
//��Ӵ�ӡ���ݵ�excal��xls,xlsx��ʽ���ݶ�ʧ��csv��ʽ���棩
//************************************
//v1.2
//����ѹ���Ǽ���
//************************************
//ͷ�ļ�����
#include "stdio.h"
#include "math.h"
//�����궨��
#define deg 360
#define pi 180
#define Pi 3.1415926
//��������
double e, r0, rt;
int Eta, Delta;
int i;
double Phi, Phis, Phi1, Phis1, dPhi;
double s0, s1, s2, s30, s31, s4, h;
double rm, Alp = 30;
double x[deg / 5 + 1], y[deg / 5 + 1], X[deg / 5 + 1], Y[deg / 5 + 1];
double A[deg / 5 + 1], a[deg / 5 + 1], pred[deg / 5 + 1];
//��������
void config();
double functions();
void theocoo(double x[], double y[]);
void pracoo(double X[], double Y[]);
double xderivation(double s, double dPhi);
double yderivation(double s, double dPhi);
void sderivation();
void printparameter();
void excalwrite();
double maxchoose();
void pressdegree();
//������
void main()
{
	config();//��ʼ��
	theocoo(x, y);//�����������
	pracoo(X, Y);//ʵ���������
	sderivation();//s����Phi��
	pressdegree();//ѹ���Ǽ���
	printparameter();//��ӡ����
	excalwrite();//д��excal
	printf_s(">�س��Լ���\n");
	getchar();//�ȴ�
}
//��ʼ������
void config()
{
	Eta = 1;//˳ʱ��
	Delta = 1;//ƫy+
	i = 0;
	dPhi = deg / 72;//72ϸ�֣�����Ϊ5
	h = 20;
	Phi = 150;
	Phis = 30;
	Phi1 = 120;
	Phis1 = 60;
	r0 = 40;
	rt = 10;
	e = 10;
	s0 = sqrt(pow(r0, 2) - pow(e, 2));
}
//s-dPhi�ֶκ���
double functions(int m, int i)
{
	//0--Phi
	s1 = (h / 2) * (1 - cos((pi / Phi) * dPhi * i / pi * Pi));
	//Phi--(Phi + Phis)
	s2 = h;
	//(Phi + Phis)--((Phi + Phis)+(Phi1/2))180-240, 36-48
	s30 = h - ((2 * h) / pow(Phi1, 2)) * pow((dPhi * i - Phi - Phis), 2);
	//((Phi + Phis + Phi1/2)--(Phi + Phis + Phi1)240-300, 48-60
	s31 = ((2 * h) / pow(Phi1, 2)) * pow((Phi + Phis + Phi1 - dPhi * i), 2);
	//(Phi + Phis + Phi1)--0
	s4 = 0;
	switch (m)
	{
	case 1: return s1; break;
	case 2: return s2; break;
	case 3: return s30; break;
	case 4: return s31; break;
	case 5: return s4; break;
	default: break;

	}
}
//����������㺯��
void theocoo(double x[], double y[])
{
	for (i = 0; i < Phi/5; i++)
	{
		x[i] = (functions(1, i) + s0) * cos(Eta * dPhi * i / pi * Pi) - Delta * e * sin(Eta * dPhi * i / pi * Pi);
		y[i] = (functions(1, i) + s0) * sin(Eta * dPhi * i / pi * Pi) + Delta * e * cos(Eta * dPhi * i / pi * Pi);
	}
	for (i; i < (Phi + Phis) / 5; i++)
	{
		x[i] = (functions(2, i) + s0) * cos(Eta * dPhi * i / pi * Pi) - Delta * e * sin(Eta * dPhi * i / pi * Pi);
		y[i] = (functions(2, i) + s0) * sin(Eta * dPhi * i / pi * Pi) + Delta * e * cos(Eta * dPhi * i / pi * Pi);
	}
	for (i; i < (Phi + Phis + Phi1 / 2) / 5; i++)
	{
		x[i] = (functions(3, i) + s0) * cos(Eta * dPhi * i / pi * Pi) - Delta * e * sin(Eta * dPhi * i / pi * Pi);
		y[i] = (functions(3, i) + s0) * sin(Eta * dPhi * i / pi * Pi) + Delta * e * cos(Eta * dPhi * i / pi * Pi);
	}
	for (i; i < (Phi + Phis + Phi1) / 5; i++)
	{
		x[i] = (functions(4, i) + s0) * cos(Eta * dPhi * i / pi * Pi) - Delta * e * sin(Eta * dPhi * i / pi * Pi);
		y[i] = (functions(4, i) + s0) * sin(Eta * dPhi * i / pi * Pi) + Delta * e * cos(Eta * dPhi * i / pi * Pi);
	}
	for (i; i <= deg / 5; i++)
	{
		x[i] = (functions(5, i) + s0) * cos(Eta * dPhi * i / pi * Pi) - Delta * e * sin(Eta * dPhi * i / pi * Pi);
		y[i] = (functions(5, i) + s0) * sin(Eta * dPhi * i / pi * Pi) + Delta * e * cos(Eta * dPhi * i / pi * Pi);
	}
	i = 0;
}
//ʵ��������㺯��(�ڰ�����)
void pracoo(double X[], double Y[])
{
	for (i = 0; i < Phi / 5; i++)
	{
		X[i] = x[i] - rt * ((yderivation(functions(1, i), dPhi * i)) / (sqrt(pow(xderivation(functions(1, i), dPhi * i), 2) + pow(yderivation(functions(1, i), dPhi * i), 2))));
		Y[i] = y[i] + rt * ((xderivation(functions(1, i), dPhi * i)) / (sqrt(pow(xderivation(functions(1, i), dPhi * i), 2) + pow(yderivation(functions(1, i), dPhi * i), 2))));
	}
	for (i; i < (Phi + Phis) / 5; i++)
	{
		X[i] = x[i] - rt * ((yderivation(functions(2, i), dPhi * i)) / (sqrt(pow(xderivation(functions(2, i), dPhi * i), 2) + pow(yderivation(functions(2, i), dPhi * i), 2))));
		Y[i] = y[i] + rt * ((xderivation(functions(2, i), dPhi * i)) / (sqrt(pow(xderivation(functions(2, i), dPhi * i), 2) + pow(yderivation(functions(2, i), dPhi * i), 2))));
	}
	for (i; i < (Phi + Phis + Phi1 / 2) / 5; i++)
	{
		X[i] = x[i] - rt * ((yderivation(functions(3, i), dPhi * i)) / (sqrt(pow(xderivation(functions(3, i), dPhi * i), 2) + pow(yderivation(functions(3, i), dPhi * i), 2))));
		Y[i] = y[i] + rt * ((xderivation(functions(3, i), dPhi * i)) / (sqrt(pow(xderivation(functions(3, i), dPhi * i), 2) + pow(yderivation(functions(3, i), dPhi * i), 2))));
	}
	for (i; i < (Phi + Phis + Phi1) / 5; i++)
	{
		X[i] = x[i] - rt * ((yderivation(functions(4, i), dPhi * i)) / (sqrt(pow(xderivation(functions(4, i), dPhi * i), 2) + pow(yderivation(functions(4, i), dPhi * i), 2))));
		Y[i] = y[i] + rt * ((xderivation(functions(4, i), dPhi * i)) / (sqrt(pow(xderivation(functions(4, i), dPhi * i), 2) + pow(yderivation(functions(4, i), dPhi * i), 2))));
	}
	for (i; i <= deg / 5; i++)
	{
		X[i] = x[i] - rt * ((yderivation(functions(5, i), dPhi * i)) / (sqrt(pow(xderivation(functions(5, i), dPhi * i), 2) + pow(yderivation(functions(5, i), dPhi * i), 2))));
		Y[i] = y[i] + rt * ((xderivation(functions(5, i), dPhi * i)) / (sqrt(pow(xderivation(functions(5, i), dPhi * i), 2) + pow(yderivation(functions(5, i), dPhi * i), 2))));
	}
	i = 0;
}
//x����dPhi��
double xderivation(double s, double dPhi)
{
	double dx, F1, F2, dd1, dd2, ddPhi;
	//x = (s + s0) * cos(Eta * dPhi) - Delta * e * sin(Eta * dPhi);
	//y = (s + s0) * sin(Eta * dPhi) + Delta * e * cos(Eta * dPhi);
	dx = 0.1;
	do{
		ddPhi = dPhi + dx;
		F1 = (s + s0) * cos(Eta * dPhi / pi * Pi) - Delta * e * sin(Eta * dPhi / pi * Pi);
		F2 = (s + s0) * cos(Eta * ddPhi / pi * Pi) - Delta * e * sin(Eta * ddPhi / pi * Pi);
		dd1 = (F2 - F1) / (dx / pi * Pi);
		dx = 0.5 * dx;
		ddPhi = dPhi + dx;
		F2 = (s + s0) * cos(Eta * ddPhi / pi * Pi) - Delta * e * sin(Eta * ddPhi / pi * Pi);
		dd2 = (F2 - F1) / (dx / pi * Pi);
	} while (fabs(dd1 - dd2 >= 1e-3));
	return dd2;
}
//y����dPhi��
double yderivation(double s, double dPhi)
{
	double dx, F1, F2, dd1, dd2, ddPhi;
	//x = (s + s0) * cos(Eta * dPhi) - Delta * e * sin(Eta * dPhi);
	//y = (s + s0) * sin(Eta * dPhi) + Delta * e * cos(Eta * dPhi);
	dx = 0.1;
	do {
		ddPhi = dPhi + dx;
		F1 = (s + s0) * sin(Eta * dPhi / pi * Pi) + Delta * e * cos(Eta * dPhi / pi * Pi);
		F2 = (s + s0) * sin(Eta * ddPhi / pi * Pi) + Delta * e * cos(Eta * ddPhi / pi * Pi);
		dd1 = (F2 - F1) / (dx / pi * Pi);
		dx = 0.5 * dx;
		ddPhi = dPhi + dx;
		F2 = (s + s0) * sin(Eta * ddPhi / pi * Pi) + Delta * e * cos(Eta * ddPhi / pi * Pi);
		dd2 = (F2 - F1) / (dx / pi * Pi);
	} while (fabs(dd1 - dd2 >= 1e-3));
	return dd2;
}
//s����dPhi��
void sderivation()
{
	int i;
	double dx, F1, F2, dd1, dd2, ddPhi;
	//0--Phi
	//s1 = (h / 2) * (1 - cos((pi / Phi) * dPhi * i / pi * Pi));
	//Phi--(Phi + Phis)
	//s2 = h;
	//(Phi + Phis)--((Phi + Phis)+(Phi1/2))180-240, 36-48
	//s30 = h - ((2 * h) / pow(Phi1, 2)) * pow((dPhi * i - Phi - Phis), 2);
	//((Phi + Phis + Phi1/2)--(Phi + Phis + Phi1)240-300, 48-60
	//s31 = ((2 * h) / pow(Phi1, 2)) * pow((Phi + Phis + Phi1 - dPhi * i), 2);
	//(Phi + Phis + Phi1)--0
	//s4 = 0;
	for (i = 0; i < Phi / 5; i++)
	{
		dx = 0.1;
		do {
			ddPhi = dPhi * i + dx;
			F1 = (h / 2) * (1 - cos((pi / Phi) * dPhi * i / pi * Pi));
			F2 = (h / 2) * (1 - cos((pi / Phi) * ddPhi / pi * Pi));
			dd1 = (F2 - F1) / (dx / pi * Pi);
			dx = 0.5 * dx;
			ddPhi = dPhi * i + dx;
			F2 = (h / 2) * (1 - cos((pi / Phi) * ddPhi / pi * Pi));
			dd2 = (F2 - F1) / (dx / pi * Pi);
		} while (fabs(dd1 - dd2 >= 1e-3));
		A[i] = fabs(dd2);
		a[i] = fabs(dd2);
	}
	for (i; i < (Phi + Phis) / 5; i++)
	{
		dx = 0.1;
		do {
			ddPhi = dPhi * i + dx;
			F1 = h;
			F2 = h;
			dd1 = (F2 - F1) / (dx / pi * Pi);
			dx = 0.5 * dx;
			ddPhi = dPhi * i + dx;
			F2 = h;
			dd2 = (F2 - F1) / (dx / pi * Pi);
		} while (fabs(dd1 - dd2 >= 1e-3));
		A[i] = fabs(dd2);
		a[i] = fabs(dd2);
	}
	for (i; i < (Phi + Phis + Phi1 / 2) / 5; i++)
	{
		dx = 0.1;
		do {
			ddPhi = dPhi * i + dx;
			F1 = h - ((2 * h) / pow(Phi1, 2)) * pow((dPhi * i - Phi - Phis), 2);
			F2 = h - ((2 * h) / pow(Phi1, 2)) * pow((ddPhi - Phi - Phis), 2);
			dd1 = (F2 - F1) / (dx / pi * Pi);
			dx = 0.5 * dx;
			ddPhi = dPhi * i + dx;
			F2 = h - ((2 * h) / pow(Phi1, 2)) * pow((ddPhi - Phi - Phis), 2);
			dd2 = (F2 - F1) / (dx / pi * Pi);
		} while (fabs(dd1 - dd2 >= 1e-3));
		A[i] = fabs(dd2);
		a[i] = fabs(dd2);
	}
	for (i; i < (Phi + Phis + Phi1) / 5; i++)
	{
		dx = 0.1;
		do {
			ddPhi = dPhi * i + dx;
			F1 = ((2 * h) / pow(Phi1, 2)) * pow((Phi + Phis + Phi1 - dPhi * i), 2);
			F2 = ((2 * h) / pow(Phi1, 2)) * pow((Phi + Phis + Phi1 - ddPhi), 2);
			dd1 = (F2 - F1) / (dx / pi * Pi);
			dx = 0.5 * dx;
			ddPhi = dPhi * i + dx;
			F2 = ((2 * h) / pow(Phi1, 2)) * pow((Phi + Phis + Phi1 - ddPhi), 2);
			dd2 = (F2 - F1) / (dx / pi * Pi);
		} while (fabs(dd1 - dd2 >= 1e-3));
		A[i] = fabs(dd2);
		a[i] = fabs(dd2);
	}
	for (i; i <= deg / 5; i++)
	{
		dx = 0.1;
		do {
			ddPhi = dPhi * i + dx;
			F1 = 0;
			F2 = 0;
			dd1 = (F2 - F1) / (dx / pi * Pi);
			dx = 0.5 * dx;
			ddPhi = dPhi * i + dx;
			F2 = 0;
			dd2 = (F2 - F1) / (dx / pi * Pi);
		} while (fabs(dd1 - dd2 >= 1e-3));
		A[i] = fabs(dd2);
	}
}
//��ӡ����
void printparameter()
{
	int j;
	printf_s("���| ���۵�\t\t\t|ʵ�ʵ�\t\t\t|ѹ����\t\t\t|ds/dPhi|(max) = %f\n", maxchoose());
	for (j = 0; j <= deg / 5; j++)
	{
		printf_s("%d: (x%d, y%d) = (%.3f, %.3f), (X%d, Y%d) = (%.3f, %.3f), Alpha = %.3f, dPhi = %d\n", j, j, j, x[j], y[j], j, j, X[j], Y[j], pred[j], 5*j);
	}
}
//�������excal��
void excalwrite()
{
	FILE *fp = NULL;
	remove(".\\zuobiao.csv");
	fp = fopen(".\\zuobiao.csv", "w");

	int i;
	fprintf(fp, "x, y, X, Y, Alpha, dPhi\n");
	for (i = 0; i <= deg/5; i++)
	{
		fprintf(fp, "%f, %f, %f, %f, %f, %d\n", x[i], y[i], X[i], Y[i], pred[i], 5*i);
	}
	fclose(fp);
	printf_s(">������д��excal, �ļ���Ϊ'zuobiao.csv'\n");
}
//ð���㷨���д�С
double maxchoose()
{
	int i, j, t;
	for (j = 0; j <= deg / 5; j++)
	{
		for (i = 0; i <= deg / 5 - 1 - j; i++)
		{
			if (a[i] < a[i + 1])
			{
				t = a[i];
				a[i] = a[i + 1];
				a[i + 1] = t;
			}
		}
	}
	return a[0];
}
//ѹ���Ǽ���
void pressdegree()
{
	int i;
	rm = (1 / atan(Alp / pi * Pi)) * maxchoose();
	for (i = 0; i <= deg / 5; i++)
	{
		pred[i] = (atan((1 / rm) * A[i])) / Pi * pi;
	}
}
