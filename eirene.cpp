
#include "eirene.hpp"
#include "eirene_impl.cpp"

using std::tuple;
using std::vector;

namespace
{
/*
    Functions that process the inital complex data.
*/

eirene::SparseBool copySubMatrix(const eirene::SparseBool A, uint32_t firstCol, uint32_t lastCol)
{
    /*
        Returns the submatrix consisting of columns [firstCol, lastCol)

        Obviously, lastCol >= firstCol.
    */
    const size_t numVal   = A.colPtrs[lastCol] - A.colPtrs[firstCol];
    const uint32_t numCol = lastCol - firstCol;

    auto rowCopies = vector<uint32_t>(numVal, 0);
    auto colCopies = vector<uint32_t>(lastCol - firstCol + 1, 0);

    for (uint32_t colInd = 0; colInd < numCol; ++colInd)
    {
        uint32_t origInd        = firstCol + colInd;
        uint32_t elemBeforeNext = A.colPtrs[origInd + 1];
        uint32_t elemBeforeCurr = A.colPtrs[origInd];

        colCopies[colInd + 1] = colCopies[colInd] + elemBeforeNext - elemBeforeCurr;

        // pointer arithmetic
        std::copy(A.rowVals.data() + elemBeforeCurr, A.rowVals.data() + elemBeforeNext,
                  rowCopies.data() + colCopies[colInd]);
    }

    return eirene::SparseBool{rowCopies, colCopies, A.numRows};
}

auto checkFiltrationMonotone(const std::vector<eirene::SparseBool>& splitComp,
                             const std::vector<std::vector<uint32_t>>& sortSplitFilt,
                             const std::vector<uint32_t>& dimPatt)

{
    bool ret = true;

    for (uint32_t dim = 0; dim < numberOfDim(dimPatt); ++dim)
    {
        const auto& dimAdj = splitComp[dim];

        for (uint32_t cell = 0; cell < dimAdj.numCols(); ++cell)
        {
            const uint32_t cellFiltVal = sortSplitFilt[dim][cell];

            for (uint32_t faceInd = dimAdj.colPtrs[cell]; faceInd < dimAdj.colPtrs[cell + 1]; ++faceInd)
            {
                uint32_t faceCell = dimAdj.rowVals[cell];
                
                // >= because they're reverse sorted, see the splitting function
                ret &= sortSplitFilt[dim - 1][faceCell] >= cellFiltVal;
            }

            if (!ret) return ret;
        }
    }

    return ret;
}

auto splitByDim(const eirene::SparseBool& complex, const vector<uint32_t>& dimPatt, const bool check = true)
{
    /*
        Splits the information up by dimension.
        Checks that bdry operator is well graded if check = true;
    */

    // Largest dimension to compute
    const uint32_t topDim = numberOfDim(dimPatt);

    vector<eirene::SparseBool> splitComplex(topDim);

    for (uint32_t currDim = 0; currDim < topDim; ++currDim)
    {
        // The columns corresponding to cells with dimension currDim are those with indexes in [dimPatt[currDim] - 1,
        // dimPatt[currDim + 1] - 1)
        // topDim is greatest lower bound on 0 chain groups... (bdry_i is 0 for i >= topDim - 1)
        if (currDim < topDim - 1)
            splitComplex[currDim] = copySubMatrix(complex, dimPatt[currDim] - 1, dimPatt[currDim + 1] - 1);

        // transform the row indexing for each dimension to be relative, i.e. start from 0 for each
        const uint32_t smallerCells        = (currDim == 0) ? 0 : dimPatt[currDim - 1] - 1,
                       smallerOrEqualCells = (currDim == topDim) ? dimPatt[currDim - 1] : dimPatt[currDim];

        for (uint32_t& rowVal : splitComplex[currDim].rowVals)
        {
            if (check && (rowVal < smallerCells || rowVal >= smallerOrEqualCells))
                return std::variant<vector<eirene::SparseBool>, std::string>{
                    "Chain boundary does not appear to be graded of dimension 1"};

            rowVal -= smallerCells;
        }

        splitComplex[currDim].numRows = smallerOrEqualCells - smallerCells - 1;
    }

    std::cerr << "Printing the split complex:\n";
    for (uint32_t currDim = 0; currDim < topDim; ++currDim)
    {
        std::cerr << "Dimension: " << currDim << "\n" << splitComplex[currDim] << "\n";
        std::cerr << "\n";
    }

    return std::variant<vector<eirene::SparseBool>, std::string>{splitComplex};
}

tuple<vector<vector<uint32_t>>, vector<float>> sortAndSplitFilts(const vector<float>& filtVals,
                                                                 const vector<uint32_t>& dimPatt)
{
    /*
        Returns a reverse sorted list of unique filtration values, and a vector of vector of indices
        into the (unique) filtration values, where we split the cells by dimension as we did with the complex
       information.
    */

    vector<vector<uint32_t>> splitIndices(numberOfDim(dimPatt));
    vector<float> uniqFilts(filtVals);

    // reverse sort, then make unique
    std::sort(uniqFilts.rbegin(), uniqFilts.rend());
    uniqFilts.erase(std::unique(uniqFilts.begin(), uniqFilts.end()), uniqFilts.end());

    // For each dim, fill a vector with the indices of the filtration value for each cell in the
    // reverse sorted list.
    uint32_t totalInd = 0;

    for (uint32_t dim = 0; dim < numberOfDim(dimPatt); ++dim)
    {
        bool overDim      = dim + 1 >= dimPatt.size();
        uint32_t numCells = overDim ? 0 : (dimPatt[dim + 1] - dimPatt[dim]);

        splitIndices[dim].resize(numCells);

        for (uint32_t cellInd = 0; cellInd < numCells; ++cellInd)
        {
            // Not using sorted property
            splitIndices[dim][cellInd] =
                std::distance(uniqFilts.begin(), std::find(uniqFilts.begin(), uniqFilts.end(), filtVals[totalInd]));
            totalInd++;
        }
    }

    return std::make_tuple(splitIndices, uniqFilts);
}

}   // namespace

std::variant<eirene::Complex, std::string> eirene::toComplex(const eirene::SparseBool& adj,
                                                             const std::vector<float>& filts,
                                                             const std::vector<uint32_t>& dimVec, const bool check)
{
    /*
     * This is really ugly because I'm playing around with using variants to fail for errors.
     * This function just splits the info by dimension and checks that it actually defines a complex.
     */

    // if (dimPatt == [])...?

    // prefix sums + 1
    auto dimPatt = computeDimensionPattern(dimVec);

    if (filts.size() != dimPatt.back() - 1)
        return std::variant<eirene::Complex, std::string>(
            "The filtration vector should have a float value for every cell.");

    // SortedSplitFilt: vector of vector of indices, each represents the filt Val of a
    // cell - so sortedSplitFilt[i][j] represents the jth cell with dimension i.
    // The index is into uniqFilt, which is a reverse sorted list of unique filtration values.
    auto [sortSplitFilt, uniqFilt] = sortAndSplitFilts(filts, dimPatt);

    // split the complex by cell dimension: in original eirene, called 'segmenting'
    auto maybeSplitComp = splitByDim(adj, dimPatt, check);

    if (std::holds_alternative<std::string>(maybeSplitComp))
        return std::variant<eirene::Complex, std::string>(std::get<std::string>(maybeSplitComp));

    const auto& splitComp = std::get<vector<eirene::SparseBool>>(maybeSplitComp);

    if (check && !checkFiltrationMonotone(splitComp, sortSplitFilt, dimPatt))
        return std::variant<eirene::Complex, std::string>("Filtration is not monotone on the cell complex poset\n");

    return std::variant<eirene::Complex, std::string>{eirene::Complex{splitComp, dimPatt, sortSplitFilt, uniqFilt}};
}

eirene::MaybeReduced eirene::eirene(const eirene::SparseBool& chainComplex, const std::vector<float>& filtVals,
                                    const std::vector<uint32_t>& dimVec)
{
    /*
        chainComplex:  A sparse matrix, representing the CW complex. Entry [i,j] is true iff cell i is
                       codimension-1 face of cell j. Thus, a row represents all the cells of which i is a (co-dim 1)
                       face, and a column represents all the the cells which are (co-dim 1) faces of cell j.
        filtVals: Filtration values for the cells
        dimVec:   Vector whose ith entry is the number of cells with dimension i ('Euler Vector')
    */

    std::variant<eirene::Complex, std::string> comp = toComplex(chainComplex, filtVals, dimVec);

    if (!std::holds_alternative<eirene::Complex>(comp))
    {
        const auto errStr = std::get<std::string>(comp);
        std::cout << errStr << std::endl;
        return MaybeReduced(errStr);
    }

    auto resultVec = eirene::performPersistence(std::get<eirene::Complex>(comp));

    return MaybeReduced(resultVec);
    // return std::variant<EireneResult, std::string>("Not done.\n");
}

