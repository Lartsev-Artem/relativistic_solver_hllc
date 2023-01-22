#include "../prj_config.h"
#ifdef SOLVE
#include "solve_global_struct.h"
#include "solve_config.h"

flux_t::flux_t()
{
	d = 0;
	v = Vector3::Zero();
	p = 0;
}

flux_t::flux_t(const Type a, const Type b, const Type c, const Type dd, const Type e)
{
	d = a;
	v = Vector3(b, c, dd);
	p = e;
}
flux_t::flux_t(const flux_t& f)
{
	d = f.d;
	v = f.v;
	p = f.p;
}

flux_t flux_t::operator+ (const flux_t& x) { d += x.d; v += x.v; p += x.p; return *this; }
flux_t flux_t::operator+= (const flux_t& x) { d += x.d; v += x.v; p += x.p; return *this; }
flux_t flux_t::operator-= (const flux_t& x) { d -= x.d; v -= x.v; p -= x.p; return *this; }
flux_t flux_t::operator* (const Type x) { d *= x; v *= x; p *= x; return *this; }
flux_t flux_t::operator- (const flux_t& x) { d -= x.d; v -= x.v; p -= x.p; return *this; }
flux_t flux_t::operator/ (const Type x) { d /= x; v /= x; p /= x; return *this; }

// ��� �������� ��� ����� �� ������ �����
Type flux_t::operator[](const int i) const 
{
	switch (i)
	{
	case 0: return d;
	case 1: return v[0];
	case 2: return v[1];
	case 3: return v[2];
	case 4: return p;

	default:
		printf("err idx in flux_t opreator []\n");
		EXIT(1);
		break;
	}
}
Type& flux_t::operator[](const int i)
{
	switch (i)
	{
	case 0: return d;
	case 1: return v[0];
	case 2: return v[1];
	case 3: return v[2];
	case 4: return p;

	default:
		printf("err idx in flux_t opreator []\n");
		EXIT(1);
		break;
	}
}
Type flux_t::operator()(const int i)
{
	switch (i)
	{
	case 0: return d;
	case 1: return v[0];
	case 2: return v[1];
	case 3: return v[2];
	case 4: return p;

	default:
		printf("err idx in flux_t opreator ()\n");
		EXIT(1);
		break;
	}
}


geo_face_t::geo_face_t()
{
	id_l = 0;
	id_r = 0;
	n = Vector3::Zero();
	S = 0;
}

geo_cell_t::geo_cell_t()
{
	for (int i = 0; i < base; i++)
	{
		id_faces[i] = 0;
		sign_n[i] = true;
	}
	V = 0;
	center = Vector3::Zero();	
}

illum_value_t::illum_value_t(const int num_dir)
{
	illum.resize(num_dir * base, 0);
	energy = 0;
	stream = Vector3::Zero();

	//Type prev_energy;
	//Vector3 prev_stream;

	impuls = Matrix3::Zero();
	div_stream = 0;
	div_impuls= Vector3::Zero();
	
	absorp_coef = 0;
	rad_en_loose_rate = 0;

	for (int i = 0; i < base; i++)
	{
		coef_inter[i] = Vector3::Zero();
	}	
}

#endif //SOVLE