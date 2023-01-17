#pragma once

#ifndef NO_EXTRAS
#include <ostream>
#endif

#include <cassert>

#define RADIANS_TO_DEGREES(A) (A * 57.295754f)
#define DEGREES_TO_RADIANS(A) (A * 0.0174533f)

namespace SPTH {
    typedef struct A2DVector {
        A2DVector () : X(0.0f), Y(0.0f) {}
        A2DVector (float GivenX, float GivenY) : X(GivenX), Y(GivenY) {}
        union {
            struct {
                float X;
                float Y;
            };
            float VectorArray[2] {};
        };
        float& operator [] (int ComponentIndex);
        A2DVector operator + (const A2DVector& AnotherVector) const;
        A2DVector operator - (const A2DVector& AnotherVector) const;
        A2DVector operator * (const A2DVector& AnotherVector) const;
        A2DVector operator * (float Constant) const;
        bool operator == (const A2DVector& AnotherVector) const;
        bool operator != (const A2DVector& AnotherVector) const;
        float Dot (const A2DVector& AnotherVector) const;
        float Magnitude () const;
        float MagnitudeSquared () const;
        float Distance (const A2DVector& AnotherVector) const;
        void Normalize ();
        A2DVector Normalized () const;
        float Angle (const A2DVector& AnotherVector) const;
        A2DVector Project (const A2DVector& TargetVector) const;
        A2DVector Perpendicular (const A2DVector& TargetVector) const;
        A2DVector Reflection (const A2DVector& NormalizedVector) const;
    } A2DVector;

    typedef struct A3DVector {
        A3DVector () : X(0.0f), Y(0.0f), Z(0.0f) {}
        A3DVector (float GivenX, float GivenY, float GivenZ) : X(GivenX), Y(GivenY) , Z(GivenZ) {}
        union {
            struct {
                float X;
                float Y;
                float Z;
            };
            float VectorArray [3] {};
        };
        float& operator [] (int ComponentIndex);
        A3DVector operator + (const A3DVector& AnotherVector) const;
        A3DVector operator - (const A3DVector& AnotherVector) const;
        A3DVector operator * (const A3DVector& AnotherVector) const;
        A3DVector operator * (float Constant) const;
        bool operator == (const A3DVector& AnotherVector) const;
        bool operator != (const A3DVector& AnotherVector) const;
        float Dot (const A3DVector& AnotherVector) const;
        A3DVector Cross (const A3DVector& AnotherVector) const;
        float Magnitude () const;
        float MagnitudeSquared () const;
        float Distance (const A3DVector& AnotherVector) const;
        void Normalize ();
        A3DVector Normalized () const;
        float Angle (const A3DVector& AnotherVector) const;
        A3DVector Project (const A3DVector& AnotherVector) const;
        A3DVector Perpendicular (const A3DVector& AnotherVector) const;
        A3DVector Reflection (const A3DVector& NormalizedVector) const;
    } A3DVector;
}

#ifndef NO_EXTRAS
std::ostream& operator << (std::ostream& OS, const SPTH::A2DVector& VectorA);
std::ostream& operator << (std::ostream& OS, const SPTH::A3DVector& VectorA);
#endif
