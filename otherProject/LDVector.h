/*
 * LLDVector.h
 *
 *  Created on: 19 may. 2023
 *      Author: ceci
 */

#ifndef LLDVector_H_
#define LLDVector_H_

#include <iostream>

using namespace std;

class LDVector
{
    private:

        unsigned asize;
    public:
        long double* vec;
        LDVector();
        LDVector(unsigned asize);
        LDVector(long double* vec, unsigned asize);
        LDVector(const LDVector &other);
        long double& operator[](unsigned index);
        const long double& operator[](unsigned index)const;


        const LDVector operator+(const LDVector &other)const;
        LDVector& operator+=(const LDVector &other);
        LDVector& operator=(const LDVector &other);
        const LDVector operator*(long double val)const;
        virtual ~LDVector();
        long double getValue(unsigned index);

        friend ostream& operator<<(ostream& os, const LDVector& vec);
        bool areEqual(const LDVector &other, double tolerance);

};





#endif /* LLDVector_H_ */
