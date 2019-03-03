#pragma once

#include <Geometry\Triangle.h>
#include <Geometry\BoundingBox.h>
#include <algorithm>
#include <ostream>
#include <cmath>
#include <string>
#include <sstream>

///////////////////////////////
// 
//////////////////////////////
Geometry::BoundingBox get_bounding_box(Geometry::Triangle const& tri)
{
	Math::Vector3f min(tri.vertex(0)), max(tri.vertex(0));
	for (unsigned int i = 0; i < 3; ++i)
	{
		min[i] = ::std::min(min[i], tri.vertex(1)[i]);
		min[i] = ::std::min(min[i], tri.vertex(2)[i]);

		max[i] = ::std::max(max[i], tri.vertex(1)[i]);
		max[i] = ::std::max(max[i], tri.vertex(2)[i]);
	}
	return Geometry::BoundingBox(min, max);
}

double my_atan(double y, double x)
{
	if (x == 0.)
	{
		return Math::piDiv2;
	}
	double res = atan(y / x);

	if (x < 0)
	{
		res -= Math::pi;
	}

	return res;
}

Math::Vector3f spherical_coordinate(Math::Vector3f const& vec)
{
	Math::Vector3f res;

	res[0] = vec.norm();

	res[1] = acos(vec[2] / res[0]);//theta: inclination

	//maybe check if x != 0
	res[2] = my_atan(vec[1], vec[0]);//phi: azimuth

	return res;
}



std::string & operator+=(std::string & str, int i)
{
	std::stringstream sstr;
	sstr << i;
	str += sstr.str();
	return str;
}

std::string & operator+=(std::string & str, double i)
{
	std::stringstream sstr;
	sstr << i;
	str += sstr.str();
	if (i == (int)i)
	{
		str += ".0";
	}
	return str;
}

std::string percent(int current, int total)
{
	std::string res;

	double p = (current * 1000.0) / total;

	p = floor(p);
	p /= 10.0;

	res += p;

	res += '%';

	return res;
}


std::string progession_bar(int current, int total, unsigned int samples = 20)
{
	std::string res = "[";
	int n = (current * samples) / total;
	res += std::string(n, '=');
	res += std::string(samples - n, ' ');
	res += "]";
	res += ": " + percent(current, total);
	return res;
}

/*
inline bool collide(Geometry::BoundingBox const& a, Geometry::BoundingBox const& b)
{
	return a.max()[0] < b.min()[0];
}
*/

/*
std::ostream & operator<<(std::ostream & out, const Geometry::Triangle & t)
{
	out << "(" << t.vertex(0) << ", " << t.vertex(1) << ", " << t.vertex(2) << ")";
	return out;
}
*/

