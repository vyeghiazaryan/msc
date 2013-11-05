/***
 * millipede: QuadEqn.cpp
 * Added by Varduhi Yeghiazaryan, 2013.
 ***/

#include "QuadEqn.h"

#include <cassert>
#include <cfloat>
#include <cmath>
#include <vector>

#include <common/math/MathConstants.h>

namespace mp {

namespace QuadEqn {

typedef std::vector<double> DoubleVector;

DoubleVector roots(double A, double B, double C)
{
	DoubleVector sol;

	assert(A > MathConstants::SMALL_EPSILON || A < -MathConstants::SMALL_EPSILON);

	double det = B * B - 4 * A * C;
	if(det > MathConstants::SMALL_EPSILON)
	{
		sol.push_back((-B - std::sqrt(det)) / (2 * A));
		sol.push_back((-B + std::sqrt(det)) / (2 * A));
	}
	else if(det > -MathConstants::SMALL_EPSILON)
	{
		sol.push_back(-B / (2 * A));
	}

	return sol;
}

double largest_root(double A, double B, double C)
{
	double largest = -DBL_MAX;
	DoubleVector sol = roots(A, B, C);
	
	for(DoubleVector::const_iterator it=sol.begin(), iend=sol.end(); it!=iend; ++it)
	{
		largest = std::max(largest, *it);
	}

	return largest;
}
double smallest_root(double A, double B, double C)
{
	double smallest = DBL_MAX;
	DoubleVector sol = roots(A, B, C);

	for(DoubleVector::const_iterator it=sol.begin(), iend=sol.end(); it!=iend; ++it)
	{
		smallest = std::min(smallest, *it);
	}
	
	return smallest;
}

}

}


