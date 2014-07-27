#ifndef _AcerDet_Core_Histogram_
#define _AcerDet_Core_Histogram_

#pragma once

#include <string>
using namespace std;

namespace AcerDet {
	namespace core {
		
		/*
		 * Container for data statistics
		 * Groups given data in separate ranges, counting
		 * amount of entries mapped to specific subrange.
		 * 
		 * Set of ranges is constructed from given meta range
		 * and number of requested subranges
		 * 
		 * Size of subrange = size of meta range / subranges count
		 * 
		 * There are two extra subranges (-inf, minimum) and (maximum, +inf)
		 * When user is inserting only porper values these two are empty.
		 */
		class Histogram {
		private:
			string title;
			double minimumValue, maximumValue;
			int resolution, storedAll, storedProperly;
			
			int* arr;
		
			inline int computePosition( double value ) const;
		public:
		
			/*
			 * Create new histogram instance
			 * 
			 * params:
			 * name - name of histogram displayed as a title when printing
			 * minVal - left edge of meta range (minimum proper value in histogram)
			 * maxVal - right edge of meta range (maximum proper value in histogram)
			 * res - histogram resolution (number of desired subranges)
			 * 
			 * throws an exception when:
			 * minVal > maxVal (inversed range is forbiden)
			 * res < 1 (at least one subrange is required)
			 */
			Histogram( const string& name, double minVal, double maxVal, unsigned int res );
			
			/*
			 * Destroy histogram and clean up all used resources
			 */ 
			~Histogram();
			
			/*
			 * Get histogram title set in constructor
			 */
			const string& getTitle() const;
			
			/*
			 * Get histogram minimum value set in constructor
			 */
			double getMinimumValue() const;
			
			/*
			 * Get histogram maximum value set in constructor
			 */
			double getMaximumValue() const;
			
			/*
			 * Get number of items registered in histogram till this call
			 */
			int storedCount() const;
			
			/*
			 * Get number of items registered properly 
			 * (between minimum and maximum)
			 * in histogram till this call
			 */
			int storedProperlyCount() const;
			
			/*
			 * Register new item
			 * value should be in range [minimumValue, maximumValue]
			 * otherwise it would be registered in one of extra subrange
			 */
			void insert( double value );
			
			/*
			 * Reset all statistics stored in histogram
			 */
			void reset();
			
			/*
			 * Print statistics to standard output (stdout by default)
			 */
			void print( bool onlyNonZero ) const;
		};
		
	}
}

#endif 
