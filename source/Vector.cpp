#include <Vector.h>

#include <cmath>
#include <cfloat>

#define COMPARISON(A, B) (fabsf((A) - (B)) <= FLT_EPSILON * fmaxf(1.0f, fmaxf(fabsf(A), fabsf(B))))

float& SPTH::A2DVector::operator [] (int ComponentIndex) {
    assert (ComponentIndex >= 0 && ComponentIndex < 2);
    return VectorArray[ComponentIndex];
}

SPTH::A2DVector SPTH::A2DVector::operator + (const A2DVector& AnotherVector) const {
    return {this->X + AnotherVector.X, this->Y + AnotherVector.Y};
}

SPTH::A2DVector SPTH::A2DVector::operator - (const A2DVector& AnotherVector) const {
    return {this->X - AnotherVector.X, this->Y - AnotherVector.Y};
}

SPTH::A2DVector SPTH::A2DVector::operator * (const A2DVector& AnotherVector) const {
    return {this->X * AnotherVector.X, this->Y * AnotherVector.Y};
}

SPTH::A2DVector SPTH::A2DVector::operator * (float Constant) const {
    return {this->X  * Constant, this->Y * Constant};
}

bool SPTH::A2DVector::operator == (const A2DVector& AnotherVector) const {
    return COMPARISON(this->X, AnotherVector.X) && COMPARISON(this->Y, AnotherVector.Y);
}

bool SPTH::A2DVector::operator != (const A2DVector& AnotherVector) const {
    return !(*this == AnotherVector);
}

float SPTH::A2DVector::Dot (const A2DVector& AnotherVector) const {
    return this->X * AnotherVector.X + this->Y * AnotherVector.Y;
}

float SPTH::A2DVector::Magnitude () const {
    return sqrtf(this->Dot(*this));
}

float SPTH::A2DVector::MagnitudeSquared () const {
    return this->Dot(*this);
}

float SPTH::A2DVector::Distance (const A2DVector& AnotherVector) const {
    A2DVector ResultVector = *this - AnotherVector;
    return ResultVector.Magnitude();
}

void SPTH::A2DVector::Normalize () {
    float Magnitude = this->Magnitude();
    assert(Magnitude != 0);
    *this = *this * (1.0f / Magnitude);
}

SPTH::A2DVector SPTH::A2DVector::Normalized () const {
    float Magnitude = this->Magnitude();
    assert(Magnitude != 0);
    return *this * (1.0f / Magnitude);
}

float SPTH::A2DVector::Angle (const A2DVector& AnotherVector) const {
    float Magnitude = sqrtf(this->MagnitudeSquared() * AnotherVector.MagnitudeSquared());
    assert(Magnitude != 0);
    return acos(this->Dot(AnotherVector) / Magnitude);
}

SPTH::A2DVector SPTH::A2DVector::Project (const A2DVector& TargetVector) const {
    float Dot = this->Dot(TargetVector);
    float MagnitudeSquared = TargetVector.MagnitudeSquared();
    assert(MagnitudeSquared != 0);
    return TargetVector * (Dot / MagnitudeSquared);
}

SPTH::A2DVector SPTH::A2DVector::Perpendicular(const SPTH::A2DVector& TargetVector) const {
    return *this - this->Project(TargetVector);
}

SPTH::A2DVector SPTH::A2DVector::Reflection (const A2DVector& NormalizedVector) const {
    float Dot = this->Dot(NormalizedVector);
    return *this - NormalizedVector * (Dot * 2.0f);
}

float& SPTH::A3DVector::operator [] (int ComponentIndex) {
    assert (ComponentIndex >= 0 && ComponentIndex < 3);
    return VectorArray[ComponentIndex];
}

SPTH::A3DVector SPTH::A3DVector::operator + (const A3DVector& AnotherVector) const {
    return {this->X + AnotherVector.X, this->Y + AnotherVector.Y, this->Z + AnotherVector.Z};
}

SPTH::A3DVector SPTH::A3DVector::operator - (const SPTH::A3DVector& AnotherVector) const {
    return {this->X - AnotherVector.X, this->Y - AnotherVector.Y, this->Z - AnotherVector.Z};
}

SPTH::A3DVector SPTH::A3DVector::operator * (const SPTH::A3DVector& AnotherVector) const {
    return {this->X * AnotherVector.X, this->Y * AnotherVector.Y, this->Z * AnotherVector.Z};
}

SPTH::A3DVector SPTH::A3DVector::operator * (float Constant) const {
    return {this->X * Constant, this->Y * Constant, this->Z * Constant};
}

bool SPTH::A3DVector::operator == (const SPTH::A3DVector& AnotherVector) const {
    return COMPARISON(this->X, AnotherVector.X) && COMPARISON(this->Y, AnotherVector.Y) && COMPARISON(this->Z, AnotherVector.Z);
}

bool SPTH::A3DVector::operator != (const SPTH::A3DVector& AnotherVector) const {
    return !(*this == AnotherVector);
}

float SPTH::A3DVector::Dot (const SPTH::A3DVector& AnotherVector) const {
    return this->X * AnotherVector.X + this->Y * AnotherVector.Y + this->Z * AnotherVector.Z;
}

SPTH::A3DVector SPTH::A3DVector::Cross (const A3DVector& AnotherVector) const {
    A3DVector ResultVector;
    ResultVector.X = this->Y * AnotherVector.Z - this->Z * AnotherVector.Y;
    ResultVector.Y = this->Z * AnotherVector.X - this->X * AnotherVector.Z;
    ResultVector.Z = this->X * AnotherVector.Y - this->Y * AnotherVector.X;
    return ResultVector;
}

float SPTH::A3DVector::Magnitude () const {
    return sqrtf(this->Dot(*this));
}

float SPTH::A3DVector::MagnitudeSquared () const {
    return this->Dot(*this);
}

float SPTH::A3DVector::Distance (const SPTH::A3DVector& AnotherVector) const {
    A3DVector ResultVector = *this - AnotherVector;
    return ResultVector.Magnitude();
}

void SPTH::A3DVector::Normalize () {
    float Magnitude = this->Magnitude();
    assert(Magnitude != 0);
    *this = *this * (1.0f / Magnitude);
}

SPTH::A3DVector SPTH::A3DVector::Normalized () const {
    float Magnitude = this->Magnitude();
    assert(Magnitude != 0);
    return *this * (1.0f / Magnitude);
}

float SPTH::A3DVector::Angle (const SPTH::A3DVector& AnotherVector) const {
    float Magnitude = sqrtf(this->MagnitudeSquared() * AnotherVector.MagnitudeSquared());
    assert(Magnitude != 0);
    return acos(this->Dot(AnotherVector) / Magnitude);
}

SPTH::A3DVector SPTH::A3DVector::Project (const SPTH::A3DVector& AnotherVector) const {
    float Dot = this->Dot(AnotherVector);
    float MagnitudeSquared = AnotherVector.MagnitudeSquared();
    assert(MagnitudeSquared != 0);
    return AnotherVector * (Dot / MagnitudeSquared);
}

SPTH::A3DVector SPTH::A3DVector::Perpendicular (const SPTH::A3DVector& AnotherVector) const {
    return *this - this->Project(AnotherVector);
}

SPTH::A3DVector SPTH::A3DVector::Reflection (const SPTH::A3DVector& NormalizedVector) const {
    float Dot = this->Dot(NormalizedVector);
    return *this - NormalizedVector * (Dot * 2.0f);
}

#ifndef NO_EXTRAS
std::ostream& operator << (std::ostream& OS, const SPTH::A2DVector& VectorA) {
    OS << "(" << VectorA.X << ", " << VectorA.Y << ")";
    return OS;
}

std::ostream& operator << (std::ostream& OS, const SPTH::A3DVector& VectorA) {
    OS << "(" << VectorA.X << ", " << VectorA.Y << ", " << VectorA.Z << ")";
    return OS;
}
#endif