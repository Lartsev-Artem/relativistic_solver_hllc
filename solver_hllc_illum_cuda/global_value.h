#ifndef GLOBAL_VALUE
#define GLOBAL_VALUE

#include <string>
#include "prj_config.h"
#include "global_def.h"

extern std::string BASE_ADRESS;


#define PI 3.14159265358979323846

const double gamma1 = 5. / 3; // 5. / 3;  // ���������� ��������
const double gamma_g = gamma1 / (gamma1 - 1);

const double eps = 1e-10;

#ifdef Sphere
const Vector3 center_point(0, 0, 0);
const Type inner_radius = 0.51; // ������ ���������� ����� (� �������)
#else 
const Vector3 center_point(10, 0, 0);
const Type inner_radius = 0.12; // ������ ���������� ����� (� �������)
#endif

//��
#define EarthMass (5.9722*1e25)
#define SunMass	  (1.9891*1e31)

// �
#define DistSun (149.6 * 10e9)
#define DistMoon 400000000.

//�/c
#define C_LIGHT 299792458.0
#define C_LIGHT_INV (1.0/(C_LIGHT))

#define DIST  (1*1e10) //DistMoon
#define MASS (1 * 1e21)//EarthMass
#define VELOCITY (3 * 1e8)//C_LIGHT

#define TIME (DIST/VELOCITY)
//#define DENSITY (MASS/(DIST*DIST*DIST))
//#define PRESSURE (MASS/(DIST*TIME*TIME))
//#define RADIATION (MASS/(TIME*TIME*TIME))
#define DENSITY (3.34*10e-24)
#define PRESSURE (DENSITY*VELOCITY*VELOCITY)
#define RADIATION (DENSITY*VELOCITY*VELOCITY*VELOCITY)

const Type R_gas = 8.314;  //������� ���������� [ ��/(����*�)]
const Type c_light = C_LIGHT; // 3 * 1e8;  //[�/c]
const Type h_plank = 6.62 * 1e-34;  // [�� * �^2 /�]
const Type k_boltzmann = 1.38 * 1e-23; // [��/K] = [ ��*�^2/(�^2*T)]
const Type sigma_thomson = 6.652 * 1e-29; //������� ������������� ��������� [m^2]
const Type m_hydrogen = 1.6735575 * 1e-27; //����� ��������[��]

#endif //GLOBAL_VALUE
