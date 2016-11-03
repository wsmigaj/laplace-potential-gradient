%module laplace_gradient
%{
#define SWIG_FILE_WITH_INIT
#include "laplace_3d_single_layer_potential_operator_gradient.hpp"
%}

%include "bempp.swg"

namespace Bempp
{
%feature("compactdefaultargs")
    laplace3dSingleLayerPotentialOperatorGradient;
} // namespace Bempp

%inline %{
namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
boost::shared_ptr<PotentialOperator<BasisFunctionType, ResultType> >
laplace3dSingleLayerPotentialOperatorGradient()
{
    typedef Bempp::Laplace3dSingleLayerPotentialOperatorGradient<BasisFunctionType, ResultType> Type;
    return boost::shared_ptr<Type>(new Type);
}

} // namespace Bempp
%}

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(laplace3dSingleLayerPotentialOperatorGradient);
}

%pythoncode %{

def createLaplace3dSingleLayerPotentialOperatorGradient(context):
    """ 
    Create and return a single-layer potential operator gradient for the Laplace
    equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object used to control the evaluation of integrals
            occurring in the definition of the potential operator.

    *Returns* a newly constructed PotentialOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.

    Note about BEM++ terminology: a *potential operator* acts on functions
    defined on a surface S and produces functions defined at any point of the
    space surrounding S, but not necessarily on S itself. In contrast, a
    *boundary operator* acts on on functions defined on a surface S and produces
    functions defined on the same surface S.
    """
    import lib
    import sys
    currentModule = sys.modules[__name__]
    basisFunctionType = context.basisFunctionType()
    resultType = context.resultType()
    result = lib._constructObjectTemplatedOnBasisAndResult(
        _laplace_gradient, "laplace3dSingleLayerPotentialOperatorGradient",
        basisFunctionType, resultType)
    result._context = context
    return result

%}
