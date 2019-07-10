// Include the (generated) file that includes the necessary data
// structures
#include "kernel_ref.h"
#include "kernel_datastructure_header.h"

#include <cassert>
#include <cmath>
#include <cstdio>

using namespace gmx;
using namespace Nbnxm;

// Help check for correct execution
template <typename T>
bool checkVector(const std::vector<T> &lhs,
                 const std::vector<T> &rhs,
                 const char *description)
{
    bool success = true;
    assert(lhs.size() == rhs.size());
    for (size_t i = 0; i != lhs.size(); ++i)
    {
        if (1e-6F < std::abs(lhs[i] - rhs[i]))
        {
            fprintf(stderr, "Failed to match %s element %zu: %g %g\n",
                    description,
                    i, lhs[i], rhs[i]);
            success = false;
        }
    }
    return success;
}

int main (int /*unused*/, char ** /*unused*/)
{
    #include "datadeclarations.h"

    // Checking for valid input
    assert(nbat.params_.comb_rule == 0);

    auto const *pairlist = &pairlistSets_.localSet_.cpuLists_[0];

    // This basic block is worth timing
    {
        nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_F_ref(pairlist,
                                                  &nbat,
                                                  &ic_,
                                                  as_rvec_array(nbat.shift_vec.data()),
                                                  &nbat.out[0]);
    }

    // Check for correct execution
    bool success = true;
    success = success && checkVector(nbat.out[0].f, out[0].f, "forces per particle and dimension");
    success = success && checkVector(nbat.out[0].fshift, out[0].fshift, "shift forces per shift vector and dimension");
    success = success && checkVector(nbat.out[0].Vvdw, out[0].Vvdw, "VDW potential per energy group");
    success = success && checkVector(nbat.out[0].Vc, out[0].Vc, "Coulomb potential per energy group");
    success = success && checkVector(nbat.out[0].VSvdw, out[0].VSvdw, "SIMD VDW potential per energy group");
    success = success && checkVector(nbat.out[0].VSc, out[0].VSc, "SIMD Coulomb potential per energy group");

    // Arrange for terminal and UNIX-style error-code feedback
    int errorCode;
    if (success)
    {
        fprintf(stdout, "Success!\n");
        errorCode = 0;
    }
    else
    {
        fprintf(stdout, "Failed!\n");
        errorCode = 1;
    }
    return errorCode;
}
