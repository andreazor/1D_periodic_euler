#include <vector> // Para std::vector
#include <omp.h>  // Para OpenMP
#include "RHSoperator.h"

template<class T>
RHSOperator<T>::RHSOperator()
{

}

template<class T>
RHSOperator<T>::~RHSOperator()
{

}

template<class T>
Central1D<T>::Central1D(DataStruct<T> &_U, DataStruct<T> &_mesh, FluxFunction<T> &_F)
    : U(_U), mesh(_mesh), F(_F)
{
    RHS.setSize(_U.getSize());
}

template<class T>
Central1D<T>::~Central1D()
{

}

template<class T>
void Central1D<T>::evalRHS(DataStruct<T> &Uin)
{
    // The BC should be included in the mesh
    // Momentarily done here by hand
    T *dataRHS = RHS.getData();
    const T *dataU = Uin.getData();
    const T *dataMesh = mesh.getData();
    const int len = U.getSize();

    // Precalcular valores de flux para evitar recalculaciones
    std::vector<T> fluxValues(len);
    for (int i = 0; i < len; ++i)
    {
        fluxValues[i] = F.computeFlux(dataU[i]);
    }

    #pragma omp parallel for
    for (int j = 1; j < len - 1; ++j)
    {
        T dx = dataMesh[j + 1] - dataMesh[j - 1];
        dataRHS[j] = -(fluxValues[j + 1] - fluxValues[j - 1]) / dx;
    }

    // Calcular el primer y œltimo elemento fuera del bucle paralelizado
    T dx0 = (dataMesh[len - 1] - dataMesh[len - 2]) + (dataMesh[1] - dataMesh[0]);
    dataRHS[0] = -(fluxValues[1] - fluxValues[len - 2]) / dx0;
    dataRHS[len - 1] = dataRHS[0]; // Aplicar condici—n de contorno
}

template<class T>
void Central1D<T>::eval()
{
    evalRHS(U);
}

template<class T>
void Central1D<T>::eval(DataStruct<T> &Uin)
{
    evalRHS(Uin);
}

template<class T>
DataStruct<T>& Central1D<T>::ref2RHS()
{
    return RHS;
}

template class RHSOperator<float>;
template class RHSOperator<double>;

template class Central1D<float>;
template class Central1D<double>;
