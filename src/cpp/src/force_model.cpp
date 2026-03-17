/**
 * HPOP Combined Force Model — Implementation
 *
 * This file provides the ForceModel base class implementations.
 * TwoBodyForce and J2Force are fully inline in force_model.h.
 *
 * The combined force model is handled by the Propagator class (propagator.cpp)
 * which iterates over a vector of ForceModel pointers and sums accelerations.
 *
 * This translation unit exists to ensure the vtable for ForceModel and its
 * concrete subclasses (DragForce, SRPForce, ThirdBodyForce) are emitted in
 * a single TU, avoiding ODR violations in multi-TU builds.
 *
 * References:
 *   - Vallado, "Fundamentals of Astrodynamics and Applications" (2013), Ch. 8
 *   - Montenbruck & Gill, "Satellite Orbits" (2000), Ch. 3
 */

#include "hpop/force_model.h"
#include "hpop/state.h"

namespace hpop {

// ForceModel destructor anchor (vtable emitted here)
// The destructor is defaulted in the header, but this TU ensures
// the vtable has a home. This is standard practice for polymorphic
// base classes with pure virtual methods.
//
// TwoBodyForce::acceleration() and J2Force::acceleration() are
// implemented inline in force_model.h for performance — they are
// called millions of times per propagation and benefit from inlining.
//
// Additional force models (DragForce, SRPForce, ThirdBodyForce) are
// implemented in their respective .cpp files:
//   - atmosphere.cpp (DragForce)
//   - srp.cpp (SRPForce)
//   - thirdbody.cpp (ThirdBodyForce)

}  // namespace hpop
