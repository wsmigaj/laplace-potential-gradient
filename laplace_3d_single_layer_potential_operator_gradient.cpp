// Copyright (C) 2016 by Wojciech Smigaj.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "laplace_3d_single_layer_potential_operator_gradient.hpp"
#include "assembly/laplace_3d_potential_operator_base_imp.hpp"

#include "fiber/explicit_instantiation.hpp"

#include "laplace_3d_single_layer_potential_gradient_kernel_functor.hpp"
#include "fiber/scalar_function_value_functor.hpp"
#include "vector_kernel_scalar_trial_integrand_functor.hpp"

#include "fiber/default_collection_of_kernels.hpp"
#include "fiber/default_collection_of_basis_transformations.hpp"
#include "fiber/default_kernel_trial_integral.hpp"

namespace Bempp
{

/** \cond PRIVATE */
template <typename BasisFunctionType, typename ResultType>
struct Laplace3dSingleLayerPotentialOperatorGradientImpl
{
    typedef Laplace3dSingleLayerPotentialOperatorGradientImpl<BasisFunctionType, ResultType>
    This;
    typedef Laplace3dPotentialOperatorBase<This, BasisFunctionType, ResultType> PotentialOperatorBase;
    typedef typename PotentialOperatorBase::KernelType KernelType;
    typedef typename PotentialOperatorBase::CoordinateType CoordinateType;

    typedef Fiber::Laplace3dSingleLayerPotentialGradientKernelFunctor<KernelType>
    KernelFunctor;
    typedef Fiber::ScalarFunctionValueFunctor<CoordinateType>
    TransformationFunctor;
    typedef Fiber::VectorKernelScalarTrialIntegrandFunctor<
    BasisFunctionType, KernelType, ResultType> IntegrandFunctor;

    Laplace3dSingleLayerPotentialOperatorGradientImpl() :
        kernels(KernelFunctor()),
        transformations(TransformationFunctor()),
        integral(IntegrandFunctor())
    {}

    Fiber::DefaultCollectionOfKernels<KernelFunctor> kernels;
    Fiber::DefaultCollectionOfBasisTransformations<TransformationFunctor>
    transformations;
    Fiber::DefaultKernelTrialIntegral<IntegrandFunctor> integral;
};
/** \endcond */

template <typename BasisFunctionType, typename ResultType>
Laplace3dSingleLayerPotentialOperatorGradient<BasisFunctionType, ResultType>::
Laplace3dSingleLayerPotentialOperatorGradient()
{
}

template <typename BasisFunctionType, typename ResultType>
Laplace3dSingleLayerPotentialOperatorGradient<BasisFunctionType, ResultType>::
~Laplace3dSingleLayerPotentialOperatorGradient()
{
}

#define INSTANTIATE_BASE_LAPLACE_SINGLE_POTENTIAL(BASIS,RESULT) \
    template class Laplace3dPotentialOperatorBase< \
Laplace3dSingleLayerPotentialOperatorGradientImpl< BASIS , RESULT >, BASIS , RESULT >
FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_BASE_LAPLACE_SINGLE_POTENTIAL);

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Laplace3dSingleLayerPotentialOperatorGradient);

} // namespace Bempp
