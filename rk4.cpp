#include "rk4.h"
#include <cassert>

template<class T>
RungeKutta4<T>::RungeKutta4(DataStruct<T> &_Un) :
Un(_Un)
{
    nSteps = 4;
    currentStep = 0;

    coeffsA = new T[4];
    coeffsB = new T[4];
    coeffsA[0] = 0.;
    coeffsA[1] = 0.5;
    coeffsA[2] = 0.5;
    coeffsA[3] = 1.;
    coeffsB[0] = 1.;
    coeffsB[1] = 2.;
    coeffsB[2] = 2.;
    coeffsB[3] = 1.;

    fi_current.setSize(Un.getSize());
    fi_accum.setSize(Un.getSize());

    Ui.setSize(Un.getSize());
}

template<class T>
RungeKutta4<T>::~RungeKutta4()
{
    delete[] coeffsA;
    delete[] coeffsB;
}

template<class T>
int RungeKutta4<T>::getNumSteps()
{
    return nSteps;
}

template<class T>
void RungeKutta4<T>::initRK()
{
    currentStep = 0;
}

template<class T>
void RungeKutta4<T>::stepUi(T dt)
{
    assert(currentStep < nSteps);

    if (currentStep == 0)
    {
        T *dataUi = Ui.getData();
        const T *dataU = Un.getData();

        for (int n = 0; n < Ui.getSize(); n++)
        {
            dataUi[n] = dataU[n];
        }
    }
    else
    {
        T *datafi = fi_current.getData();
        T *dataUi = Ui.getData();
        const T *dataU = Un.getData();

        for (int n = 0; n < Ui.getSize(); n++)
        {
            dataUi[n] = dataU[n] + coeffsA[currentStep] * dt * datafi[n];
        }
    }
}

template<class T>
void RungeKutta4<T>::finalizeRK(const T dt)
{
    T *dataUn = Un.getData();
    T *dataUi = Ui.getData();
    T *dataFi_accum = fi_accum.getData();

    // set Ui to 0
    for (int n = 0; n < Ui.getSize(); n++)
    {
        dataUi[n] = dataFi_accum[n];
    }

    const T oneDiv6 = 1. / 6.;
    for (int n = 0; n < Ui.getSize(); n++)
    {
        dataUn[n] += dt * oneDiv6 * dataUi[n];
    }
}

template<class T>
void RungeKutta4<T>::setFi(DataStruct<T> &_F)
{
    T *dataFi_current = fi_current.getData();
    T *dataFi_accum = fi_accum.getData();
    const T *dataF = _F.getData();
    const T b = coeffsB[currentStep];

    for (int n = 0; n < Ui.getSize(); n++)
    {
        dataFi_current[n] = dataF[n];
        dataFi_accum[n] += b * dataF[n];
    }

    currentStep++;
}

template<class T>
DataStruct<T> *RungeKutta4<T>::currentU()
{
    return &Ui;
}

template class RungeKutta4<float>;
template class RungeKutta4<double>;
