#include "Histogram.h"
#include <cstdio>
#include <cmath>

using namespace AcerDet::core; 

Histogram::Histogram( const string& name, double minVal, double maxVal, unsigned int res ) :
	title( name ),
	minimumValue( minVal ),
	maximumValue( maxVal ),
	resolution( res ),
	storedAll(0),
	storedProperly(0)
{
	// wrong interval
	if (maxVal < minVal) {
		printf ("Acer DET 2.0 Histogram ctor -> end of given interval is less than its begin!\n");
		throw -1;
	}
	
	// zero resoulution
	if (res < 1) {
		printf ("Acer DET 2.0 Histogram ctor -> cannot create a histogram without resolution!\n");
		throw -1;
	}
	
	arr = new int[res + 2];
	for (unsigned int i=0;i<res+2;++i)
		arr[i] = 0;
} 

Histogram::~Histogram() {
	delete[] arr;
}

/*
 * [0] - less than minimum
 * [1,...,res] - normal range
 * [res+1] - greater than maximum
 */
#define epsi 0.000001
int Histogram::computePosition( double val ) const {
	if (val < minimumValue)
		return 0;
	
	if (val > maximumValue)
		return resolution+1;
		
	double coef = (val - minimumValue) / (maximumValue - minimumValue);
	if (coef > 1.0 - epsi) return resolution;
	return 1 + floor(coef * (double)resolution);
}

const string& Histogram::getTitle() const { return title; }
double Histogram::getMinimumValue() const { return minimumValue; }
double Histogram::getMaximumValue() const { return maximumValue; }
int Histogram::storedCount() const { return storedAll; }
int Histogram::storedProperlyCount() const { return storedProperly; }

void Histogram::insert( double value ) {
	int insertionIndex = computePosition(value);
	arr[ insertionIndex ]++;
	
	if (insertionIndex != 0 && insertionIndex != resolution+1)
		storedProperly++;
	
	storedAll++;
}

void Histogram::reset() {
	for (int i=0;i<resolution+2;++i)
		arr[i] = 0;
		
	storedAll = storedProperly = 0;
}

void Histogram::print( bool onlyNonZero ) const {
	printf ("AcerDet Histogram: \"%s\"\n", title.c_str());
	
	if (arr[0] > 0) 
		printf ("(-inf, %lf) = %d\n", minimumValue, arr[0]);
	
	for (int i=0;i<resolution;++i) {
		if (arr[i+1] == 0)
			continue;
			
		double coef_lower = (double)i / (double)resolution;
		double coef_upper = (double)(i+1) / (double)resolution;
		double lower = minimumValue + (maximumValue - minimumValue) * coef_lower;
		double upper = minimumValue + (maximumValue - minimumValue) * coef_upper;
		printf ("[%lf, %lf) = %d\n", lower, upper, arr[1+i]);
	}
	
	if (arr[resolution+1] > 0)
		printf ("(%lf, +inf) = %d\n", maximumValue, arr[resolution+1]);
}
