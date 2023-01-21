#ifndef GLOBAL_VALUE
#define GLOBAL_VALUE

#include <string>
#include "prj_config.h"
#include "global_def.h"

extern std::string BASE_ADRESS;


#define PI 3.14159265358979323846

const double gamma1 = 5. / 3; // 5. / 3;  // показатель адиабаты
const double gamma_g = gamma1 / (gamma1 - 1);

const double eps = 1e-10;

#define C_LIGHT 299792458.0
#define C_LIGHT_INV (1.0/(C_LIGHT))

#ifdef Sphere
const Vector3 center_point(0, 0, 0);
const Type inner_radius = 0.51; // радиус внутренней сферы (с запасом)
#else 
const Vector3 center_point(10, 0, 0);
const Type inner_radius = 0.12; // радиус внутренней сферы (с запасом)
#endif

#endif //GLOBAL_VALUE
