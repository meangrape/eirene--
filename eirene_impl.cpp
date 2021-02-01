#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <iostream>

using std::tuple;
using std::vector;

namespace
{
auto reindex(const std::vector<uint32_t>& indexVec, const std::vector<uint32_t>& sourceVec)
{
    uint32_t outSz = std::min(indexVec.size(), sourceVec.size());

    vector<uint32_t> result(outSz);

    for (uint32_t i = 0; i < outSz; ++i) result[i] = sourceVec[indexVec[i]];

    return result;
}

inline uint32_t numberOfDim(const vector<uint32_t>& dimPatt)
{
    /*
     * If the complex has maximum cell dimension of n, returns n + 1.
     */
    if (dimPatt.size() == 0) { std::cerr << "Error! dimPatt should never have size 0!" << std::endl; }

    return dimPatt.size();
}

inline int32_t quickSearchSorted(const vector<uint32_t>& vec, uint32_t key, uint32_t start, uint32_t stop)
{
    for (uint32_t i = start; (vec[i] <= key) && (i < stop); ++i)
        if (vec[i] == key) return i;

    return -1;
}

void applyPermutation(vector<uint32_t>& vals, vector<uint32_t>& inds)
{
    /*
        Applies the permutation specified by inds to vals, in-place.
        Overwrites both.
    */
    for (uint32_t i = 0; i < inds.size(); ++i)
    {
        uint32_t holeInd = i;

        while (i != inds[holeInd])
        {
            uint32_t next = inds[holeInd];

            std::swap(vals[holeInd], vals[next]);
            inds[holeInd] = holeInd;
            holeInd       = next;
        }

        inds[holeInd] = holeInd;
    }
}

vector<uint32_t> getSortPerm(const vector<float>& fltVec, const bool rev = false)
{
    /*
        Returns the permutation that sortes fltVec, if rev then for reverse order.
    */
    auto temp = vector<tuple<float, uint32_t>>(fltVec.size());

    for (uint32_t i = 0; i < fltVec.size(); ++i) temp[i] = std::make_tuple(fltVec[i], i);

    if (rev)
        std::sort(temp.rbegin(), temp.rend());   // reverse sort
    else
        std::sort(temp.begin(), temp.end());

    auto toRet = vector<uint32_t>(fltVec.size());
    std::transform(temp.begin(), temp.end(), toRet.begin(), [](auto tup) { return std::get<uint32_t>(tup); });

    return toRet;
}

vector<uint32_t> uniqueIntsNotInUnsortedUpTo(const vector<uint32_t>& uniqIntVec, uint32_t endPt,
                                             const uint32_t subVecEnd = 0)
{
    /*
        Returns the integers from the interval [0, endPt) NOT present in
        uniqIntVec[0 : subVecEnd>0 ? subVecEnd : uniqIntVec.size()].

        Note that as suggested, uniqIntVec may be unsorted vector of unique nonnnegative integers.

        This only works if endPt >= uniqIntVec.size(), and (x \in uniqIntVec) => (x < endPt)
    */

    const uint32_t vecSz = subVecEnd > 0 ? subVecEnd : uniqIntVec.size();
    vector<uint32_t> complement{};

    // If vecSz is 0, act like every integer in [0, endPt) is not present,
    // so the complement is the whole interval
    if (vecSz == 0)
    {
        complement.resize(endPt);
        std::iota(complement.begin(), complement.end(), 0);
        return complement;
    }

    // Else, if vecSz == endPt, act like evry integer is present, so
    // the complement is nothing
    else if (vecSz == endPt)
    {
        return complement;
    }

    else if (endPt < vecSz)
    {
        // will error below
        std::cerr << "Error! endPt < size of passed vector in uniqueIntsNotInUnsortedUpTo!\n" << std::endl;
    }

    vector<bool> complementSupport(endPt, true);
    for (uint32_t i = 0; i < vecSz; ++i) complementSupport[uniqIntVec[i]] = false;

    complement.resize(endPt - vecSz);
    uint32_t marker = 0;

    for (uint32_t i = 0; i < endPt; ++i)
    {
        if (complementSupport[i])
        {
            complement[marker] = i;
            marker++;
        }
    }

    if (marker != (endPt - vecSz))
    {
        std::cerr << "Stop here\n";
        std::cerr << "Warning! uniqueIntsNotInUnsortedUpTo is not working as expected..." << std::endl;
        std::getchar();
    }

    return complement;
}

inline vector<uint32_t> getColFiltTransitionPtr(const vector<uint32_t>& colFilt, const uint32_t upTo)
{
    /*
        ColFilt is a sorted list of integers, this returns a vector whose ith element is the
        beginning of the ith string of constant values in colFilt.
    */

    vector<uint32_t> toRet(0);

    if (colFilt.size() == 0) return toRet;

    toRet.reserve(colFilt.size() / 2);   // random guess
    toRet.push_back(0);
    uint32_t currVal = colFilt[0];

    for (uint32_t i = 1; i < std::min(static_cast<size_t>(upTo), colFilt.size()); ++i)
    {
        if (colFilt[i] != currVal)
        {
            toRet.push_back(i);
            currVal = colFilt[i];
        }
    }

    toRet.push_back(upTo + 1);
    return toRet;
}

vector<uint32_t> someWeirdPerm(const vector<uint32_t>& coBdryColPtrs, const vector<uint32_t>& columnFiltPtr,
                               uint32_t numLowCells)
{
    /*
        Given that ColumnFiltPtr is a sorted list of indices
        Returns a permutation vector such that:
           1) For every i < numLowCells, the range columnFiltPtr[i]:columnFiltPtr[i + 1] maps to itself under the
           permutation
           2) For every i < numLowCells, the indices in columnFiltPtr[i]:columnFiltPtr[i] correspond to columns
           in coBdryPtrs in an order such that the number of elements in the indexed columns is nondecreasing.
    */

    // doing std::sort() on each column's indices with the sorting function being coBdryColPtrs[i + 1]
    // -coBdryColPtrs[i]. Also just sorting the whole thing ignoring criticalLowCount
    // make it fast sometime

    vector<uint32_t> toRet(numLowCells);
    std::iota(toRet.begin(), toRet.end(), 0);

    const auto sortFunc = [&coBdryColPtrs](uint32_t colIndA, uint32_t colIndB) {
        uint32_t aNumElem = coBdryColPtrs[colIndA + 1] - coBdryColPtrs[colIndA];
        uint32_t bNumElem = coBdryColPtrs[colIndB + 1] - coBdryColPtrs[colIndB];

        return aNumElem < bNumElem;
    };

    for (uint32_t colFiltValInd = 0; colFiltValInd < columnFiltPtr.size(); ++colFiltValInd)
    {
        uint32_t endDist = std::min(columnFiltPtr[colFiltValInd + 1], numLowCells);
        auto start = toRet.begin() + columnFiltPtr[colFiltValInd], end = toRet.begin() + endDist;
        std::sort(start, end, sortFunc);

        if (endDist == numLowCells) break;
    }

    return toRet;
}

void getPairs(eirene::SparseBool& coBdry, const vector<uint32_t>& rowFilt, const vector<uint32_t>& colFilt,
              const uint32_t criticalHighCount, const uint32_t criticalLowCount, eirene::PairVec& pairVec)
{
    /*
        Finds pairs for the filt/morse function.
        As usual, coBdry is highCell x lowCell
    */

    vector<uint32_t>&coBdryRowVals = coBdry.rowVals, &coBdryColPtrs = coBdry.colPtrs;
    vector<uint32_t> largestFiltCoFace(criticalLowCount, 0);
    vector<uint32_t> rowWiseSum(criticalHighCount, 0);
    pairVec.clear();

    /*
        For each low dim. cell (coBdry col), if it has cofaces, find the
        one with the largest filtration index (corresp. to low birth time of coface?). Also as we go through the
        cofaces (rows), aggregate the sum of all values in each row (number of cells with that row's cell as a coface)
    */

    for (uint32_t lowCellInd = 0; lowCellInd < criticalLowCount; ++lowCellInd)
    {
        uint32_t highestCoFaceInd = coBdryColPtrs[lowCellInd],
                 highestFiltVal   = rowFilt[coBdryRowVals[highestCoFaceInd]];

        // look at cofaces
        for (uint32_t coFaceInd = coBdryColPtrs[lowCellInd] + 1; coFaceInd < coBdryColPtrs[lowCellInd + 1]; ++coFaceInd)
        {
            uint32_t highCellInd  = coBdryRowVals[coFaceInd];
            uint32_t otherFiltVal = rowFilt[highCellInd];

            ++rowWiseSum[highCellInd];

            if (otherFiltVal > highestFiltVal)
            {
                highestCoFaceInd = coFaceInd;
                highestFiltVal   = otherFiltVal;
            }
        }

        largestFiltCoFace[lowCellInd] = highestCoFaceInd;
    }

    vector<uint32_t> columnFiltPtr        = getColFiltTransitionPtr(colFilt, criticalLowCount);
    vector<uint32_t> colWiseSumLinearized = someWeirdPerm(coBdryColPtrs, columnFiltPtr, criticalLowCount);
    vector<uint32_t> colsInNameOrder(criticalLowCount);
    std::iota(colsInNameOrder.begin(), colsInNameOrder.end(), 0);
    applyPermutation(colsInNameOrder, colWiseSumLinearized);

    vector<bool> nCoveredSupp(criticalHighCount, true);

    // Are cols sorted such that first indices after perm have small coface degree??

    for (uint32_t colPermIndex = 0; colPermIndex < criticalLowCount; ++colPermIndex)
    {
        const uint32_t colIndex   = colsInNameOrder[colPermIndex];
        const uint32_t numCoFaces = coBdryColPtrs[colIndex + 1] - coBdryColPtrs[colIndex];

        if (numCoFaces > 0)
        {
            uint32_t largeFiltFaceInd = largestFiltCoFace[colIndex],       // index in sparse bdry
                largeCoFace           = coBdryRowVals[largeFiltFaceInd],   // row index of coface
                largeCoFaceDegree     = rowWiseSum[largeCoFace],           // num low cells with coface as coface
                largeCoFaceFiltVal    = rowFilt[largeCoFace];

            for (uint32_t otherCoFaceInd = largeFiltFaceInd + 1; otherCoFaceInd < coBdryColPtrs[colIndex + 1];
                 ++otherCoFaceInd)
            {
                uint32_t otherCoFace = coBdryRowVals[otherCoFaceInd], otherFiltVal = rowFilt[otherCoFace],
                         otherDegree = rowWiseSum[otherCoFace];

                if (otherFiltVal == largeCoFaceFiltVal && otherDegree <= largeCoFaceDegree &&
                    !nCoveredSupp[largeCoFace] && nCoveredSupp[otherCoFace])
                {
                    largeCoFaceDegree = otherDegree;
                    largeCoFace       = otherCoFace;
                }
            }

            if (nCoveredSupp[largeCoFace]) pairVec.push_back(eirene::MorsePair{colIndex, largeCoFace});

            for (uint32_t colElem = coBdryColPtrs[colIndex]; colElem < coBdryColPtrs[colIndex + 1]; ++colElem)
            {
                nCoveredSupp[coBdryRowVals[colElem]] = false;
            }
        }
    }
}

}   // anonymous namespace

// sparse linear algebra
namespace spla
{
struct F2Result
{
    /*
     * Vectors of this struct are used in sparse F2 matrix operations.
     * Allows merging three length "numRows" arrays into one.
     * Makes using it a little more confusing but simplifies cache locality and
     * reduces # of vectors created.
     */
    uint32_t colInd;
    uint32_t rowInd;
    bool indicator;
};

inline void F2MultAlternative(const uint32_t rowInd, const uint32_t colInd, uint32_t& prodRowCounter,
                              vector<F2Result>& colResult)
{
    /*
     *  Forms the computation of the inner loop of many sparse F2 algorithms.
     */
    if (colResult[rowInd].colInd != (colInd + 1))
    {
        colResult[rowInd].colInd         = colInd + 1;
        colResult[rowInd].indicator      = true;
        colResult[prodRowCounter].rowInd = rowInd;
        prodRowCounter++;
    }
    else
        colResult[rowInd].indicator = !colResult[rowInd].indicator;
}

inline void F2SilentPremult(const uint32_t startInd, const uint32_t endInd, const uint32_t colInd,
                            const vector<uint32_t>& rowVals, vector<F2Result>& colResult)
{
    /*
     *  A "Silent" sparse matrix has implicit 1s on the diagonal.
     *  When doing operations on these matrices, this gets called in the beginning of
     *  each loop computing a column of the output matrix.
     */

    for (uint32_t valInd = startInd; valInd < endInd; valInd++)
    {
        uint32_t rowInd = rowVals[valInd];

        colResult[valInd - startInd].rowInd = rowInd;
        colResult[rowInd].colInd            = colInd + 1;
        colResult[rowInd].indicator         = true;
    }
}

auto F2InverseSilentOut(const eirene::SparseBool& A)
{
    /*
     * Returns some sort of silent inverse of a matrix.
     */

    const uint32_t aNumCols = A.numCols(), aNumRows = A.numRows;

    eirene::SparseBool aInv{aNumCols, aNumRows};

    aInv.rowVals.reserve(A.rowVals.size());
    vector<F2Result> colResult(aNumCols);

    for (uint32_t aColInd = 0; aColInd < aNumCols; ++aColInd)
    {
        // If the column has only one nonzero element, the inverse has none
        if (A.colPtrs[aColInd] + 1 == A.colPtrs[aColInd + 1]) aInv.colPtrs[aColInd + 1] = aInv.colPtrs[aColInd];

        // Elif the columns has two, so that the high dimension cell has two cofaces...
        // All sparse should be sorted now, need to test this...
        else if (A.colPtrs[aColInd] + 2 == A.colPtrs[aColInd + 1])
        {
            uint32_t k = 0;

            // if the row index of the element is < column index,
            // set k to that row index. Otherwise, set k to the
            // row index of the next one.
            // This is the unsorted rows bit. Note we're taking as assumption that
            // one of the elements in the columns satisfies row index < column index,
            // makes sense from a boundary perspective
            if (A.rowVals[A.colPtrs[aColInd]] < aColInd)
                k = A.rowVals[A.colPtrs[aColInd]];
            else
                k = A.rowVals[A.colPtrs[aColInd] + 1];

            // append row indices of elements in column k
            aInv.rowVals.insert(aInv.rowVals.end(), aInv.rowVals.begin() + aInv.colPtrs[k],
                                aInv.rowVals.begin() + aInv.colPtrs[k + 1]);

            aInv.rowVals.push_back(k);
            aInv.colPtrs[aColInd + 1] = aInv.rowVals.size();
        }

        // General case:
        else
        {
            F2SilentPremult(A.colPtrs[aColInd], A.colPtrs[aColInd + 1], aColInd, A.rowVals, colResult);

            colResult[aColInd].indicator = false;

            uint32_t newRowCounter = A.colPtrs[aColInd + 1] - A.colPtrs[aColInd];

            for (uint32_t colElemInd = A.colPtrs[aColInd]; colElemInd < A.colPtrs[aColInd + 1]; ++colElemInd)
            {
                uint32_t elemRowVal = A.rowVals[colElemInd];

                // unsorted rows again
                if (elemRowVal < aColInd)
                {
                    for (uint32_t aInvColElInd = aInv.colPtrs[elemRowVal]; aInvColElInd < aInv.colPtrs[elemRowVal + 1];
                         ++aInvColElInd)
                    {
                        uint32_t k = aInv.rowVals[aInvColElInd];

                        F2MultAlternative(k, aColInd, newRowCounter, colResult);
                    }
                }
            }

            aInv.rowVals.reserve(aInv.rowVals.size() + newRowCounter);

            for (uint32_t newRowInd = 0; newRowInd < newRowCounter; ++newRowInd)
            {
                uint32_t rowInd = colResult[newRowInd].rowInd;

                if (colResult[rowInd].indicator) aInv.rowVals.push_back(rowInd);
            }

            aInv.colPtrs[aColInd + 1] = aInv.rowVals.size();
            std::sort(aInv.rowVals.begin() + aInv.colPtrs[aColInd], aInv.rowVals.begin() + aInv.colPtrs[aColInd + 1]);
        }
    }

    return aInv;
}

auto F2InverseSilentInAndOut(const eirene::SparseBool& A)
{
    /*
        Returns the "morse inverse" of a matrix.
        Has something to do with only caring about entries i,j with i < j (upper triangular)
    */

    const uint32_t aNumCols = A.numCols();

    // It seems a morse inverse assumes that the input matrix is square?
    // eirene::SparseBool aInv{aNumCols, aNumRows};
    eirene::SparseBool aInv{aNumCols, aNumCols};

    aInv.rowVals.reserve(A.rowVals.size());
    vector<F2Result> colResult(aNumCols);

    for (uint32_t aColInd = 0; aColInd < aNumCols; ++aColInd)
    {
        // If the column has no one nonzero elements, the inverse has none
        if (A.colPtrs[aColInd] == A.colPtrs[aColInd + 1]) aInv.colPtrs[aColInd + 1] = aInv.colPtrs[aColInd];

        // Elif the columns has one
        else if (A.colPtrs[aColInd] + 1 == A.colPtrs[aColInd + 1])
        {
            uint32_t k = A.rowVals[A.colPtrs[aColInd]];

            // append indices of elements in column k
            aInv.rowVals.insert(aInv.rowVals.end(), aInv.rowVals.begin() + aInv.colPtrs[k],
                                aInv.rowVals.begin() + aInv.colPtrs[k + 1]);

            aInv.rowVals.push_back(k);
            aInv.colPtrs[aColInd + 1] = aInv.rowVals.size();
        }

        // General case:
        else
        {
            F2SilentPremult(A.colPtrs[aColInd], A.colPtrs[aColInd + 1], aColInd, A.rowVals, colResult);

            uint32_t newRowCounter = A.colPtrs[aColInd + 1] - A.colPtrs[aColInd];

            for (uint32_t colElemInd = A.colPtrs[aColInd]; colElemInd < A.colPtrs[aColInd + 1]; ++colElemInd)
            {
                uint32_t elemRowVal = A.rowVals[colElemInd];

                for (uint32_t aInvColElInd = aInv.colPtrs[elemRowVal]; aInvColElInd < aInv.colPtrs[elemRowVal + 1];
                     ++aInvColElInd)
                {
                    uint32_t k = aInv.rowVals[aInvColElInd];

                    F2MultAlternative(k, aColInd, newRowCounter, colResult);
                }
            }

            aInv.rowVals.reserve(aInv.rowVals.size() + newRowCounter);

            for (uint32_t newRowInd = 0; newRowInd < newRowCounter; ++newRowInd)
            {
                uint32_t rowInd = colResult[newRowInd].rowInd;

                if (colResult[rowInd].indicator) aInv.rowVals.push_back(rowInd);
            }

            aInv.colPtrs[aColInd + 1] = aInv.rowVals.size();
        }
    }

    return aInv;
}

auto matMultF2(const eirene::SparseBool& A, const eirene::SparseBool& B)
{
    const size_t aNumRows = A.numRows, bNumCols = B.numCols();

    eirene::SparseBool prodMat{aNumRows, bNumCols};

    size_t sizeGuess = std::min(aNumRows * bNumCols, static_cast<size_t>(A.numNonZero() + B.numNonZero()));
    prodMat.rowVals.reserve(sizeGuess);

    vector<F2Result> colResult(aNumRows);

    for (uint32_t prodColInd = 0; prodColInd < bNumCols; ++prodColInd)
    {
        uint32_t newRowCounter = 0;

        for (uint32_t bValInd = B.colPtrs[prodColInd]; bValInd < B.colPtrs[prodColInd + 1]; ++bValInd)
        {
            uint32_t bRowInd = B.rowVals[bValInd];

            for (uint32_t aValInd = A.colPtrs[bRowInd]; aValInd < A.colPtrs[bRowInd + 1]; ++aValInd)
            {
                uint32_t aRowInd = A.rowVals[aValInd];

                F2MultAlternative(aRowInd, prodColInd, newRowCounter, colResult);
            }
        }

        /*        uint32_t prodColRows = std::accumulate(colResult.begin(), colResult.end(), 0,
                                                       [](uint32_t acc, const F2Result& res) { return acc +
           res.indicator; });

                const uint32_t totalNonZero     = prodMat.colPtrs[prodColInd] + prodColRows;
                prodMat.colPtrs[prodColInd + 1] = totalNonZero;

                if (prodColRows == 0) continue;

                prodMat.rowVals.reserve(totalNonZero);

                for (uint32_t j = 0; j < prodColRows; ++j)
                {
                    uint32_t rowVal = colResult[j].rowInd;

                    if (colResult[rowVal].indicator) prodMat.rowVals.push_back(rowVal);
                }

                std::sort(prodMat.rowVals.begin() + prodMat.colPtrs[prodColInd], prodMat.rowVals.begin() +
           totalNonZero);
          */

        prodMat.colPtrs[prodColInd + 1] = prodMat.colPtrs[prodColInd];

        if (newRowCounter == 0) continue;

        for (uint32_t j = 0; j < newRowCounter; ++j)
        {
            uint32_t rowVal = colResult[j].rowInd;

            if (colResult[rowVal].indicator) prodMat.rowVals.push_back(rowVal);
        }

        prodMat.colPtrs[prodColInd + 1] = prodMat.rowVals.size();
    }

    return prodMat;
}

auto silentLeftMatMultF2(const eirene::SparseBool& A, const eirene::SparseBool& B)
{
    /*
        A * B = C over F2 for sparse matrices, with A "silent".

        Silent means: implicit 1s on the diagonal.

        Output prodMAt is not silent.

        c_ij = \Sigma_k a_ik * b_kj
    */

    const uint32_t aNumRows = A.numRows, bNumCols = B.numCols();

    eirene::SparseBool prodMat{aNumRows, bNumCols};
    prodMat.rowVals.reserve(A.colPtrs.back() + B.colPtrs.back());

    vector<F2Result> colResult(aNumRows);

    for (uint32_t bColInd = 0; bColInd < bNumCols; ++bColInd)
    {
        // Each iteration computes a column of the product,
        uint32_t prodNonzeroCounter = B.colPtrs[bColInd + 1] - B.colPtrs[bColInd];

        // "Silent in"
        F2SilentPremult(B.colPtrs[bColInd], B.colPtrs[bColInd + 1], bColInd, B.rowVals, colResult);

        for (uint32_t bValInd = B.colPtrs[bColInd]; bValInd < B.colPtrs[bColInd + 1]; ++bValInd)
        {
            uint32_t bRowInd = B.rowVals[bValInd];

            for (uint32_t aValInd = A.colPtrs[bRowInd]; aValInd < A.colPtrs[bRowInd + 1]; aValInd++)
            {
                // A has nonzero entry at (aRowInd, bRowInd)
                uint32_t aRowInd = A.rowVals[aValInd];
                F2MultAlternative(aRowInd, bColInd, prodNonzeroCounter, colResult);
            }
        }

        prodMat.colPtrs[bColInd + 1] = prodMat.colPtrs[bColInd];

        if (prodNonzeroCounter == 0) continue;

        for (uint32_t j = 0; j < prodNonzeroCounter; ++j)
        {
            uint32_t rowVal = colResult[j].rowInd;

            if (colResult[rowVal].indicator) prodMat.rowVals.push_back(rowVal);
        }

        prodMat.colPtrs[bColInd + 1] = prodMat.rowVals.size();
    }

    return prodMat;
}

void testSilents()
{
    eirene::SparseBool A(3, 3), B(3, 3);

    for (uint32_t i = 0; i < 3; ++i)
    {
        if (i == 2) A.emplace_at(1, i - 1, i);
        ;

        B.emplace_at(1, i, i);
        ;
    }

    std::cerr << "silent left test: \n";
    std::cerr << "A:\n" << A << "\nB:\n" << B << "\nC:\n\n";

    auto C = spla::silentLeftMatMultF2(A, B);

    std::cerr << C << std::endl;
    std::getchar();
}

auto blockProdSum(const eirene::SparseBool& C, const eirene::SparseBool& B, const eirene::SparseBool& D,
                  eirene::SparseBool& result)
{
    /*
     *  Calculates S = D + C*B, stores the result in prod, overwrites prod!
     */

    // note ?? below
    bool validSizes = (D.numRows == C.numRows) && (C.numCols() == B.numRows) && (B.numCols() == D.numCols());
    if (!validSizes) std::cerr << "Invalid sizes in blockProdSum\n" << std::endl;

    const uint32_t dNumCols = D.numCols(), dNumRows = D.numRows;

    result.colPtrs.resize(dNumCols + 1);
    result.rowVals.clear();
    std::fill(result.colPtrs.begin(), result.colPtrs.end(), 0);
    result.numRows = D.numRows;

    vector<F2Result> colResult(dNumRows);

    for (uint32_t resultColInd = 0; resultColInd < dNumCols; ++resultColInd)
    {
        F2SilentPremult(D.colPtrs[resultColInd], D.colPtrs[resultColInd + 1], resultColInd, D.rowVals, colResult);

        uint32_t newRowsCounter = D.colPtrs[resultColInd + 1] - D.colPtrs[resultColInd];

        // ??
        if (B.colPtrs[resultColInd] < B.colPtrs[dNumCols])
        {
            for (uint32_t bElemInd = B.colPtrs[resultColInd]; bElemInd < B.colPtrs[resultColInd + 1]; ++bElemInd)
            {
                uint32_t bRowInd = B.rowVals[bElemInd];

                for (uint32_t cElemInd = C.colPtrs[bRowInd]; cElemInd < C.colPtrs[bRowInd + 1]; ++cElemInd)
                {
                    uint32_t cRowVal = C.rowVals[cElemInd];

                    F2MultAlternative(cRowVal, resultColInd, newRowsCounter, colResult);
                }
            }
        }

        result.colPtrs[resultColInd + 1] = result.colPtrs[resultColInd];

        if (newRowsCounter == 0) continue;

        for (uint32_t j = 0; j < newRowsCounter; ++j)
        {
            uint32_t rowVal = colResult[j].rowInd;

            if (colResult[rowVal].indicator) result.rowVals.push_back(rowVal);
        }

        result.colPtrs[resultColInd + 1] = result.rowVals.size();
    }
}

auto blockProdSumSilentIColsLeftWrite(const eirene::SparseBool& C, const eirene::SparseBool& B,
                                      const eirene::SparseBool& D, eirene::SparseBool& result,
                                      const vector<uint32_t>& col2SilentI)
{
    /*
     *  Calculates S = D + C*B, stores the result in prod, overwrites prod!
     *  Does something with col2SilentI tho, here prod is the transform/index tracking matrix...
     *  col2SilentI[j] = "critical cell in the ith pair column??" Doesn't make sense
     */

    const uint32_t dNumRows = D.numRows, dNumCols = D.numCols();

    result.colPtrs.resize(D.numCols() + 1);
    std::fill(result.colPtrs.begin(), result.colPtrs.end(), 0);
    result.rowVals.clear();
    result.numRows = D.numRows;

    vector<F2Result> colResult(dNumRows);

    for (uint32_t resultColInd = 0; resultColInd < dNumCols; ++resultColInd)
    {
        F2SilentPremult(D.colPtrs[resultColInd], D.colPtrs[resultColInd + 1], resultColInd, D.rowVals, colResult);

        uint32_t newRowsCounter = D.colPtrs[resultColInd + 1] - D.colPtrs[resultColInd];

        for (uint32_t bElemInd = B.colPtrs[resultColInd]; bElemInd < B.colPtrs[resultColInd + 1]; ++bElemInd)
        {
            uint32_t bRowInd = B.rowVals[bElemInd];
            uint32_t iRowVal = col2SilentI[bRowInd];

            F2MultAlternative(iRowVal, resultColInd, newRowsCounter, colResult);
        }

        for (uint32_t bElemInd = B.colPtrs[resultColInd]; bElemInd < B.colPtrs[resultColInd + 1]; ++bElemInd)
        {
            uint32_t bRowInd = B.rowVals[bElemInd];

            for (uint32_t cElemInd = C.colPtrs[bRowInd]; cElemInd < C.colPtrs[bRowInd + 1]; ++cElemInd)
            {
                uint32_t cRowVal = C.rowVals[cElemInd];

                F2MultAlternative(cRowVal, resultColInd, newRowsCounter, colResult);
            }
        }

        result.colPtrs[resultColInd + 1] = result.colPtrs[resultColInd];

        if (newRowsCounter == 0) continue;

        for (uint32_t j = 0; j < newRowsCounter; ++j)
        {
            uint32_t rowVal = colResult[j].rowInd;

            if (colResult[rowVal].indicator) result.rowVals.push_back(rowVal);
        }

        result.colPtrs[resultColInd + 1] = result.rowVals.size();
    }
}

auto stackedSubMat(const eirene::SparseBool& blockMat, const vector<uint32_t>& rows1, const vector<uint32_t>& rows2,
                   const vector<uint32_t>& cols)
{
    /*
        Returns the sparse submatrices whose elements are in columns cols, and rows rows 1/2.

        The row values in the submatrices are sorted and the actual row indices are relative to the
        order in rows1.

        Rows1 .> Rows2?
    */

    const uint32_t numCols = cols.size(), numRows = blockMat.numRows;
    const vector<uint32_t>&colPtrs = blockMat.colPtrs, &rowVals = blockMat.rowVals;

    auto rowSupp = vector<std::tuple<uint8_t, uint32_t>>(numRows, std::make_tuple(0, 0));

    for (uint32_t i = 0; i < rows1.size(); ++i) rowSupp[rows1[i]] = std::make_tuple(1, i);
    for (uint32_t i = 0; i < rows2.size(); ++i) rowSupp[rows2[i]] = std::make_tuple(2, i);

    std::array<uint32_t, 3> nonZero{0, 0, 0};

    for (uint32_t rowVal : blockMat.rowVals) ++nonZero[std::get<uint8_t>(rowSupp[rowVal])];

    // Create the submatrix rows/cols:
    eirene::SparseBool topBlock{rows1.size(), numCols}, bottomBlock{rows2.size(), numCols};
    topBlock.rowVals.reserve(nonZero[1]);
    bottomBlock.rowVals.reserve(nonZero[2]);

    for (uint32_t col = 0; col < cols.size(); ++col)   // colInd : cols)
    {
        uint32_t colInd = cols[col];

        for (uint32_t elemInd = colPtrs[colInd]; elemInd < colPtrs[colInd + 1]; ++elemInd)
        {
            uint32_t rowVal    = rowVals[elemInd];
            uint8_t subMatInd  = std::get<uint8_t>(rowSupp[rowVal]);
            uint32_t subRowInd = std::get<uint32_t>(rowSupp[rowVal]);

            if (subMatInd == 1)
            {
                auto spot = std::lower_bound(topBlock.rowVals.begin() + topBlock.colPtrs[col], topBlock.rowVals.end(),
                                             subRowInd);
                topBlock.rowVals.insert(spot, subRowInd);
            }
            else if (subMatInd == 2)
            {
                auto spot = std::lower_bound(bottomBlock.rowVals.begin() + bottomBlock.colPtrs[col],
                                             bottomBlock.rowVals.end(), subRowInd);
                bottomBlock.rowVals.insert(spot, subRowInd);
            }
        }

        topBlock.colPtrs[col + 1]    = topBlock.rowVals.size();
        bottomBlock.colPtrs[col + 1] = bottomBlock.rowVals.size();
    }

    return std::make_tuple(topBlock, bottomBlock);
}

vector<uint32_t> findDownStreamInUpperTriangMat(const eirene::SparseBool& coBdry, const vector<uint32_t>& unpairedCols,
                                                const eirene::PairVec& pairVec)
{
    /*
     * Returns ("keeps") pairs whose columns do have a 1 in the same row as currently unpaired
     * columns. This is mentioned in the literatrue, for example see the paper "Duality in Persistent Homology",
     * section 4.2 titled "Practice".
     */

    const uint32_t numPairs = pairVec.size(), numRows = coBdry.numRows;
    const auto &coBdryColPtrs = coBdry.colPtrs, &coBdryRowVals = coBdry.rowVals;

    if (numPairs == 0) return vector<uint32_t>{};

    // Find out what rows are nonzero in the columns in coBoundary with index in unpairedCols
    // This might be faster with a set or sorted array instead of two passes
    vector<bool> criticalCellRowSupp(numRows);

    for (uint32_t colInd : unpairedCols)
        for (uint32_t elemInd = coBdryColPtrs[colInd]; elemInd < coBdryColPtrs[colInd + 1]; ++elemInd)
            criticalCellRowSupp[coBdryRowVals[elemInd]] = true;

    vector<uint32_t> rowSup(0);
    rowSup.reserve(numRows);

    for (uint32_t i = 0; i < numRows; ++i)
        if (criticalCellRowSupp[i]) rowSup.push_back(i);

    // gives the pair represented by the row, given the row
    // rowTranslator[rowIndex] = pairIndex
    vector<uint32_t> rowTranslator = vector<uint32_t>(numRows);
    for (uint32_t i = 0; i < numPairs; ++i) rowTranslator[pairVec[i].highInd] = i;

    // given a row index, returns true iff a pair is supported by that row
    vector<bool> pairedRowsupp(numRows, false);
    for (uint32_t i = 0; i < numPairs; ++i) pairedRowsupp[pairVec[i].highInd] = true;

    // For each coface of a critical/unpaired lower dimensional cell,
    // mark any pairs which share that coface
    // dSS[pairInd] = true iff a critical column has a 1 in it's row
    vector<bool> downStreamSupp(pairVec.size(), false);
    for (uint32_t rowInd : rowSup)
        if (pairedRowsupp[rowInd]) downStreamSupp[rowTranslator[rowInd]] = true;

    // For each pair, backwards...
    for (int32_t pairInd = numPairs - 1; pairInd >= 0; --pairInd)
    {
        // For each coface of the low dim cell in the pair...
        const uint32_t pairCol = pairVec[pairInd].lowInd;

        for (uint32_t elemInd = coBdryColPtrs[pairCol]; elemInd < coBdryColPtrs[pairCol + 1]; ++elemInd)
        {
            uint32_t rowVal = coBdryRowVals[elemInd];

            // if a pair is supported on the row (for example, when rowVal == pairRows[pairInd])
            // AND that pair shares a coface with a critical lower dimensional cell AKA
            // there's a 1 in the row in the critical cell part of the matrix
            if (pairedRowsupp[rowVal] && downStreamSupp[rowTranslator[rowVal]])
            {
                // then go through all of the elements in the column again
                for (uint32_t elemInd = coBdryColPtrs[pairCol]; elemInd < coBdryColPtrs[pairCol + 1]; ++elemInd)
                {
                    uint32_t rowVal = coBdryRowVals[elemInd];

                    // And mark any paired cofaces as having a critical face cell
                    if (pairedRowsupp[rowVal]) downStreamSupp[rowTranslator[rowVal]] = true;
                }

                break;
            }
        }
    }

    // Finally, return the indexes of pairs which share share cofaces with critical cells
    const uint32_t downStreamNum = std::accumulate(downStreamSupp.begin(), downStreamSupp.end(), 0);

    vector<uint32_t> downStreamPairInds(downStreamNum);
    uint32_t counter = 0;

    // i iterates over pairs
    for (uint32_t i = 0; counter < downStreamNum; ++i)
        if (downStreamSupp[i]) downStreamPairInds[counter++] = i;

    return downStreamPairInds;
}

void insertCols(const eirene::SparseBool& A, eirene::SparseBool& B, const vector<uint32_t>& columnInds,
                const uint32_t destInd)
{
    /*
        Copies columns with indices in columnInds from sparse matrix A to B,
        starting at column destInd, and overwriting the existing columns in B.
    */
    uint32_t numNewElem = 0;
    for (uint32_t ind : columnInds) numNewElem += A.colPtrs[ind + 1] - A.colPtrs[ind];

    B.rowVals.resize(B.colPtrs[destInd] + numNewElem);

    if (B.colPtrs.size() < destInd + columnInds.size() + 1) B.colPtrs.resize(destInd + columnInds.size() + 1);

    B.colPtrs[0] = 0;   // ?

    for (uint32_t newCol = 0; newCol < columnInds.size(); ++newCol)
    {
        uint32_t newColInd = destInd + newCol, oldColInd = columnInds[newCol];
        B.colPtrs[newColInd + 1] = B.colPtrs[newColInd] + A.colPtrs[oldColInd + 1] - A.colPtrs[oldColInd];
        std::copy(A.rowVals.begin() + A.colPtrs[oldColInd], A.rowVals.begin() + A.colPtrs[oldColInd + 1],
                  B.rowVals.begin() + B.colPtrs[newColInd]);
    }
}

inline static auto copySubMatrix(const eirene::SparseBool& mat, const vector<uint32_t>& colIndices)
{
    /*
        Returns a pair of vectors, representing the sparse submatrix [:] x [colIndices[:]).
    */

    uint32_t numValAcc = 0;
    for (uint32_t colInd = 0; colInd < colIndices.size(); ++colInd)
        numValAcc += mat.colPtrs[colIndices[colInd] + 1] - mat.colPtrs[colIndices[colInd]];

    const uint32_t numCol = colIndices.size(), numVal = numValAcc;

    eirene::SparseBool subMat{mat.numRows, numCol};
    subMat.rowVals.resize(numVal);

    for (uint32_t colInd = 0; colInd < numCol; ++colInd)
    {
        uint32_t origInd        = colIndices[colInd];
        uint32_t elemBeforeNext = mat.colPtrs[origInd + 1];
        uint32_t elemBeforeCurr = mat.colPtrs[origInd];

        subMat.colPtrs[colInd + 1] = subMat.colPtrs[colInd] + elemBeforeNext - elemBeforeCurr;

        // pointer arithmetic
        std::copy(mat.rowVals.data() + elemBeforeCurr, mat.rowVals.data() + elemBeforeNext,
                  subMat.rowVals.data() + subMat.colPtrs[colInd]);
    }

    return subMat;
}

void schurIt(eirene::SparseBool& coBdry, vector<uint32_t>& criticalLow, vector<uint32_t>& criticalHigh,
             eirene::PairVec& pairVec, eirene::PairVec& schurPairs, const vector<uint32_t>& unpairedRows,
             const vector<uint32_t>& unpairedCols, eirene::SparseBool& transMat, eirene::SparseBool& schurForm)
{
    /*
     * Does an iteration of Block Schur decomposition.
     */

    // Record some values from before Schur
    const uint32_t numPairs = pairVec.size();

    // If there are new pairs since last iteration, need to do bookkeeping
    vector<uint32_t> indsToCopy(numPairs);
    std::transform(pairVec.rbegin(), pairVec.rend(), indsToCopy.begin(),
                   [](const eirene::MorsePair& pair) { return pair.lowInd; });

    insertCols(transMat, schurForm, indsToCopy, schurPairs.size());

    std::transform(pairVec.rbegin(), pairVec.rend(), std::back_inserter(schurPairs),
                   [&criticalHigh, &criticalLow](const eirene::MorsePair& pair) {
                       return eirene::MorsePair{criticalLow[pair.lowInd], criticalHigh[pair.highInd]};
                   });

    // Returns a vector of indices of pairs in pairVec which share cofaces with critical cells.
    // Pairs whose columns do not share cofaces can be "thrown away", reducing the size of remaining
    // matrix to reduce.
    auto pairsKept = findDownStreamInUpperTriangMat(coBdry, unpairedCols, pairVec);

    const uint32_t numKept = pairsKept.size();
    auto pairedRows = vector<uint32_t>(numKept), pairedCols = vector<uint32_t>(numKept);

    for (uint32_t i = 0; i < numKept; ++i)
    {
        uint32_t keptIndex = pairsKept[i];
        pairedRows[i]      = pairVec[keptIndex].highInd;
        pairedCols[i]      = pairVec[keptIndex].lowInd;
    }

    /*
    pairedCols unpairedCols
    _____ _______

    [ A .. - B ..]  |
    [ ..   - ..  ]  | pairedRows
    [----------- ]
    [C ..  - D ..]  |
    [..    - ..  ]  | unpairedRows

    Note pairedRows/Cols contains only those that were kept by downstream, the rest disappear.
    */

    auto [A, C] = stackedSubMat(coBdry, pairedRows, unpairedRows, pairedCols);
    auto [B, D] = stackedSubMat(coBdry, pairedRows, unpairedRows, unpairedCols);

    auto L = copySubMatrix(transMat, pairedCols);
    auto R = copySubMatrix(transMat, unpairedCols);

    // Get "silent/morse" inverse of the paired cell / top left submatrix
    // Result is pairedCols x pairedRows
    auto Ai = F2InverseSilentOut(A);

    // Multiply aInv x B, dimensions are (pairedCols x pairedRows) x (pairedRows x unpairedCols),
    // so results in E of size pairedCols x unpairedCols.
    // E not silent!
    auto E = silentLeftMatMultF2(Ai, B);

    // Translates a kept pair's index to
    vector<uint32_t> translator(numKept);
    for (uint32_t j = 0; j < numKept; ++j) translator[j] = criticalLow[pairedCols[j]];

    // Updates the row/col to high/low cell mapping
    for (uint32_t i = 0; i < unpairedRows.size(); ++i) criticalHigh[i] = criticalHigh[unpairedRows[i]];
    criticalHigh.resize(unpairedRows.size());

    for (uint32_t i = 0; i < unpairedCols.size(); ++i) criticalLow[i] = criticalLow[unpairedCols[i]];
    criticalLow.resize(unpairedCols.size());

    // at this point, low/high cell ref reflect not the dimensions of coBdry but of submatrix D (unpaired x
    // unpaired) Also, dimension: C : (unpairedRows x pairedCols), E : (pairedCols x unpairedCols), D :
    // (unpairedRows x unpairedCols)
    //  ! This is an inplace operation on coBdry !
    // Computes D + C * E
    blockProdSum(C, E, D, coBdry);

    // updates translator
    blockProdSumSilentIColsLeftWrite(L, E, R, transMat, translator);
}

auto morseLU(eirene::SparseBool& coBdry, eirene::PairVec& pairVec, const vector<uint32_t>& highFiltTemp,
             const vector<uint32_t>& lowFiltTemp)
{
    /*
        Here bdryRow/Col is a sparse representation of the coboundary.
        Thus, the rows represent the higher dimensional cells, and the cols the lower.
        The convention below is that coboundary is of size highNum x lowNum;
        In the original eirene, this was m x n;

        criticalLow contains the indices of lower dimensional cells that were unpaired in the last iteration.
        criticalHigh contains the indices of all high dimensional cells.
        Both of these indexes are defined relative to the unpaired cells passed into this algorithm, as
        opposed to the global indexes in the complex, and the cells they refer to in lowlab (not passed in)
        are in sorted order w.r.t the filtration.

        Then, unpairedRows/Cols defined below indicate the row/cols that correspond to these critical high/low
        cells. Thus it ends up partially representing the permutation that reduces the matrix.

        As the matrix is reduced, it shrinks, and so do criticalLow/Hi and unpairedRows/Cols.
    */
    vector<uint32_t>& coBdryColPtrs = coBdry.colPtrs;

    // uint32_t criticalLowCount = criticalLow.size(), criticalHighCount = criticalHigh.size(),
    uint32_t criticalLowCount = coBdry.numCols(), criticalHighCount = coBdry.numRows,
             nonZeroVals = coBdryColPtrs[criticalLowCount] + 1;

    vector<uint32_t> criticalHigh(criticalHighCount), criticalLow(criticalLowCount);
    std::iota(criticalLow.begin(), criticalLow.end(), 0);
    std::iota(criticalHigh.begin(), criticalHigh.end(), 0);

    // Set up the indexing so that the paired rows and cols come first,
    // followed by the critical cells in unpairedrows/cols
    uint32_t maxNumPairs = std::min(criticalLowCount, criticalHighCount), numPairs = pairVec.size();

    pairVec.reserve(maxNumPairs);

    auto schurPairs = eirene::PairVec();
    schurPairs.reserve(maxNumPairs);

    auto unpairedRows = vector<uint32_t>(criticalHighCount - numPairs),
         unpairedCols = vector<uint32_t>(criticalLowCount - numPairs);

    std::iota(unpairedRows.begin(), unpairedRows.end(), numPairs);
    std::iota(unpairedCols.begin(), unpairedCols.end(), numPairs);

    eirene::SparseBool transform{criticalHighCount, criticalLowCount}, schurForm{criticalHighCount, criticalLowCount};

    schurIt(coBdry, criticalLow, criticalHigh, pairVec, schurPairs, unpairedRows, unpairedCols, transform, schurForm);

    nonZeroVals = std::max(nonZeroVals, coBdryColPtrs[criticalLowCount] + 1);

    auto rowFilt = vector<uint32_t>(unpairedRows.size()), colFilt = vector<uint32_t>(unpairedCols.size());

    uint32_t counter = 0;
    while (coBdryColPtrs[criticalLowCount] > 0)
    {
        counter++;

        std::transform(criticalHigh.begin(), criticalHigh.end(), rowFilt.begin(),
                       [&highFiltTemp](const uint32_t cellInd) { return highFiltTemp[cellInd]; });
        std::transform(criticalLow.begin(), criticalLow.end(), colFilt.begin(),
                       [&lowFiltTemp](const uint32_t cellInd) { return lowFiltTemp[cellInd]; });

        // finds new pairVec
        getPairs(coBdry, rowFilt, colFilt, criticalHighCount, criticalLowCount, pairVec);

        // Just find the unpaired Rows and Columns as those that are not in the list of paired row/columns.
        // Should merge this into one algo.
        const uint32_t numPairs = pairVec.size();

        vector<uint32_t> halfPairVec(numPairs);
        std::transform(pairVec.begin(), pairVec.end(), halfPairVec.begin(),
                       [](const eirene::MorsePair pair) { return pair.highInd; });
        unpairedRows = uniqueIntsNotInUnsortedUpTo(halfPairVec, criticalHigh.size(), numPairs);

        std::transform(pairVec.begin(), pairVec.end(), halfPairVec.begin(),
                       [](const eirene::MorsePair pair) { return pair.lowInd; });
        unpairedCols      = uniqueIntsNotInUnsortedUpTo(halfPairVec, criticalLow.size(), numPairs);
        criticalLowCount  = unpairedCols.size();
        criticalHighCount = unpairedRows.size();

        // CriticalHigh and criticalLow are modified at the end of schurit, set to size unpairedCols, which is const in
        // schurIt
        schurIt(coBdry, criticalLow, criticalHigh, pairVec, schurPairs, unpairedRows, unpairedCols, transform,
                schurForm);

        nonZeroVals = std::max(nonZeroVals, coBdryColPtrs.back() + 1);
    }

    uint32_t lastSRowMarker = schurForm.colPtrs[schurPairs.size()],
             lastTRowMarker = transform.colPtrs[criticalLowCount];

    // truncate the schur form and transform
    schurForm.rowVals.resize(lastSRowMarker);
    schurForm.colPtrs.resize(schurPairs.size());
    transform.rowVals.resize(lastTRowMarker);
    transform.colPtrs.resize(criticalLowCount + 1);

    std::transform(transform.colPtrs.begin(), transform.colPtrs.end(), transform.colPtrs.begin(),
                   [lastSRowMarker](uint32_t val) { return val + lastSRowMarker; });

    // Append the transform to the Schur matrix (really more of an addition to a blank part of schurForm)
    schurForm.colPtrs.insert(schurForm.colPtrs.end(), transform.colPtrs.begin(), transform.colPtrs.end());
    schurForm.rowVals.insert(schurForm.rowVals.end(), transform.rowVals.begin(), transform.rowVals.end());

    std::vector<uint32_t> tlab(0);
    tlab.reserve(schurPairs.size());
    std::transform(schurPairs.begin(), schurPairs.end(), std::back_inserter(tlab),
                   [](const eirene::MorsePair& pair) { return pair.lowInd; });
    tlab.insert(tlab.end(), criticalLow.begin(), criticalLow.end());

    return std::make_tuple(schurForm, schurPairs, tlab);
}

auto transposeSparseSubmatrix(const eirene::SparseBool& sparseMat, const vector<uint32_t>& rows,
                              const vector<uint32_t>& cols)
{
    /*
        Returns the sparse representation of the transpose of the submatrix given delineated
        by the arguments.

        The outputs' rows are in the order given in cols. Does not assume rows, cols are sorted.

        Notation: input matrix has dimension numRows x (colPtrs.size() - 1)
                  output (transposed) submatrix has dimensions cols.size() x rows.size() = subCols x subRows
    */
    const vector<uint32_t>&rowVals = sparseMat.rowVals, &colPtrs = sparseMat.colPtrs;
    const uint32_t subRows = rows.size();

    eirene::SparseBool transMat{cols.size(), rows.size()};
    transMat.colPtrs.resize(subRows + 1);
    transMat.numRows = cols.size();

    // Compute the column counts of transposed submatrix:
    for (uint32_t elemInd = 0; elemInd < rowVals.size(); ++elemInd)
    {
        auto subRowPos = std::find(rows.begin(), rows.end(), rowVals[elemInd]);
        if (subRowPos == rows.end()) continue;

        uint32_t newRowInd = std::distance(rows.begin(), subRowPos);
        for (uint32_t j = newRowInd + 1; j < subRows + 1; ++j) transMat.colPtrs[j]++;
    }

    const uint32_t numNotZero = transMat.colPtrs.back();
    transMat.rowVals.reserve(numNotZero);

    // Sort the entries into their final resting spots (bad algo)
    // Should this not reindex the rows as well?
    for (uint32_t rowInd : rows)       // for each of the rows (transposed cols)
        for (uint32_t colInd : cols)   // go through the cols (transposed rows) in the order given
            for (uint32_t colSt = colPtrs[colInd]; colSt < colPtrs[colInd + 1]; ++colSt)
                if (rowVals[colSt] == rowInd) transMat.rowVals.push_back(colInd);
    // assert(transMat.rowVals.size() == sparseMat.rowVals.size());
    return transMat;
}

auto maxNonsingularBlock(const eirene::SparseBool& mat, const eirene::PairVec& pairVec, const vector<uint32_t>& tid,
                         const uint32_t numLowerPairs)
{
    /*
     *  Returns the submatrix of mat, which is the original boundary operator, with rows in tid and
     *  columns in highPairInds, in the order they are in in those vectors.
     *
     *  Length of tid is equal to numCritialLow.
     */
    const uint32_t numLowCells = mat.numRows, numPairs = pairVec.size(), numCriticalLow = numLowCells - numLowerPairs,
                   numRowsOut = tid.size();

    if (numPairs == 0) return eirene::SparseBool(numRowsOut, numCriticalLow);

    vector<uint32_t> emptyVec{};
    vector<uint32_t> highPairInds(numPairs);

    std::transform(pairVec.begin(), pairVec.end(), highPairInds.begin(), [](const auto& pair) { return pair.highInd; });

    auto [A, Z] = stackedSubMat(mat, tid, emptyVec, highPairInds);

    const uint32_t numElem = A.numNonZero();
    A.colPtrs.resize(numCriticalLow + 1);
    std::fill(A.colPtrs.begin() + numPairs, A.colPtrs.end(), numElem);

    return A;
}

}   // namespace spla

namespace eirene
{
auto getCycleRep(const eirene::SparseBool& lowBoundary, const std::array<SparseBool, 3>& factors,
                 const std::array<SparseBool, 3>& lowFactors, const PairVec& pairs, const vector<uint32_t>& tid,
                 const vector<uint32_t>& lowTid, const uint32_t twoBelowCells, const uint32_t cycleNumber)
{
    /*
     * Li has the information for the (co)cycle group.
     *
     *
     *
     */
    // Note it only uses Li from the highFactors
    const auto& [L, R, Li] = factors;

    const uint32_t numChainSummands = Li.colPtrs[cycleNumber + 1] - Li.colPtrs[cycleNumber];

    // Get the cycles?
    // tid reindexed row indexes from Li give cocycle cell indexes
    vector<uint32_t> summands(numChainSummands + 1);

    for (uint32_t i = 0; i < numChainSummands; ++i)
    {
        auto LiInd  = Li.colPtrs[cycleNumber] + i;
        summands[i] = tid[Li.rowVals[LiInd]];
    }
    summands.back() = tid[cycleNumber];

    // Any coboundaries?
    if (twoBelowCells == 0) return summands;

    // If so, need to compute cocycles mod coboundaires
    // Good thing we factorized the (dim - 2) -> (dim - 1) coboundary
    auto support                    = vector<bool>(twoBelowCells, false);
    const auto& [lowL, lowR, lowLi] = lowFactors;

    // pass the (dim - 1) cells in summands into the (dim - 1) -> (dim - 2) boundary
    // So support has the boundary of our cocycle Li
    for (auto j : summands)
    {
        for (uint32_t k = lowBoundary.colPtrs[j]; k < lowBoundary.colPtrs[j + 1]; ++k)
        {
            uint32_t i = lowBoundary.rowVals[k];
            support[i] = !support[i];
        }
    }

    uint32_t suppCard = std::accumulate(support.begin(), support.end(), 0);

    auto brv = vector<uint32_t>(0);
    brv.reserve(suppCard);

    // get those values from support that lie in lowTid into brv
    for (uint32_t lowTidIndex = 0; lowTidIndex < lowTid.size(); ++lowTidIndex)
    {
        uint32_t lowIndex = lowTid[lowTidIndex];
        if (support[lowIndex]) brv.push_back(lowTidIndex);
    }

    // all in one column
    auto bcp = vector<uint32_t>{0, static_cast<uint32_t>(brv.size())};
    auto b   = eirene::SparseBool{brv, bcp, L.numCols()};

    // b = R(L(b))
    // "RLM is the fundamental cycle matrix of of the d-dimensional boundary operator with respect
    // to the doubly-minimal basis we have constructed, where M the submatrix of
    // the boundary with rows indexed by tid." - GHP

    b = spla::silentLeftMatMultF2(lowL, b);
    b = spla::silentLeftMatMultF2(lowR, b);

    vector<uint32_t> lowToHigh(twoBelowCells, 0);
    // lowToHigh[pairVec[i].lowInd] = pairVec[i].highInd
    for (uint32_t i = 0; i < pairs.size(); ++i) lowToHigh[pairs[i].lowInd] = pairs[i].highInd;

    // "Nonzero entries might lie in nonbasis rows"
    for (uint32_t i = 0; i < b.rowVals.size(); ++i) b.rowVals[i] = lowToHigh[lowTid[b.rowVals[i]]];
    b.rowVals.insert(b.rowVals.end(), summands.begin(), summands.end());

    return b.rowVals;
}

auto getCycleReps(const PairVec& pairs, const PairVec& lowPairs, const vector<uint32_t>& tid,
                  const vector<uint32_t>& lowTid, const vector<std::array<eirene::SparseBool, 3>>& factors,
                  const SparseBool& lowBoundary, const vector<uint32_t>& lowFilt, const vector<uint32_t>& highFilt,
                  const uint32_t twoBelowCells, const uint32_t dim)
{
    vector<vector<uint32_t>> cycles;
    cycles.reserve(tid.size());

    // mortal cycles, from pairs
    for (uint32_t pairInd = 0; pairInd < pairs.size(); ++pairInd)
    {
        auto pair = pairs[pairInd];
        if (lowFilt[pair.lowInd] == highFilt[pair.highInd]) continue;
        
        cycles.emplace_back(
            getCycleRep(lowBoundary, factors[dim], factors[dim - 1], lowPairs, tid, lowTid, twoBelowCells, pairInd));
    }

    // immortals
    for (uint32_t cycleInd = pairs.size(); cycleInd < tid.size(); ++cycleInd)
        cycles.emplace_back(
            getCycleRep(lowBoundary, factors[dim], factors[dim - 1], lowPairs, tid, lowTid, twoBelowCells, cycleInd));

    return cycles;
}

auto extractHomology(vector<MorseReducedResult>& reducedVec, const vector<SparseBool>& splitComplex,
                     vector<uint32_t> dimPatt, const vector<vector<uint32_t>>& sortedSplitFilt)
{
    /*
     *  After doing the Morse LU transform on all the boundary matrices, calculate homology.
     *
     *  From orig eirene:
     *      L, Li, and R are all indexed by tid, just like (trv,tcp)
     *      L (not Li) is the transpose of (trv,tcp)
     *      up to permutation and deletion of some zero rows, RLM is the fundamental
     *      cycle matrix of the d-dimensional boundary operator with respect
     *      to the doubly-minimal basis we have constructed, where M the submatrix of
     *      the boundary with rows indexed by tid.
     */

    vector<std::array<eirene::SparseBool, 3>> dimFactors(numberOfDim(dimPatt));
    const uint32_t topDim = numberOfDim(dimPatt);

    for (uint32_t dim = 1; dim < numberOfDim(dimPatt); ++dim)
    {
        const uint32_t numLowCells    = splitComplex[dim - 1].numCols(),
                       numCriticalLow = numLowCells - reducedVec[dim - 1].pairVec.size();

        SparseBool reindexedReduced{reducedVec[dim].LUCoBdry};

        vector<uint32_t> lowTranslator(numLowCells, 0);
        auto& tid = reducedVec[dim].tid;

        for (uint32_t i = 0; i < tid.size(); ++i) lowTranslator[tid[i]] = i;

        // rowVals[i] = lowTranslator[rowVals[i]]
        reindexedReduced.rowVals = reindex(reindexedReduced.rowVals, lowTranslator);

        vector<uint32_t> rowsOfInterest(numCriticalLow);
        std::iota(rowsOfInterest.begin(), rowsOfInterest.end(), 0);

        // Could save a little time by first taking transpose of reindexedReduced,
        // except that F2Inverse only works on upper triangular mats
        auto Li = spla::F2InverseSilentInAndOut(reindexedReduced);

        vector<uint32_t> cols(Li.numCols());
        std::iota(cols.begin(), cols.end(), 0);

        Li = spla::transposeSparseSubmatrix(Li, rowsOfInterest, cols);

        // for the last dim, only want to calculate Li (only boundaries/cocycles)
        if (dim == topDim)
        {
            dimFactors[dim][2] = Li;
            continue;
        }

        cols.resize(reindexedReduced.numCols());
        std::iota(cols.begin(), cols.end(), 0);

        // As noted in the source, the first (rank of M) elements of tid index the
        // complete set of nonzero columns in the reduced matrix
        auto nonSingular = spla::maxNonsingularBlock(splitComplex[dim], reducedVec[dim].pairVec, reducedVec[dim].tid,
                                                     reducedVec[dim - 1].pairVec.size());

        auto L = spla::transposeSparseSubmatrix(reindexedReduced, rowsOfInterest, cols);
        auto b = spla::silentLeftMatMultF2(L, nonSingular);
        auto R = spla::F2InverseSilentOut(b);

        dimFactors[dim] = std::array<SparseBool, 3>{L, R, Li};
    }

    for (uint32_t dim = 1; dim < numberOfDim(dimPatt); ++dim)
    {
        const uint32_t twoBelowCells = dim > 1 ? splitComplex[dim - 2].numCols() : 0;

        reducedVec[dim].cycleReps = getCycleReps(
            reducedVec[dim].pairVec, reducedVec[dim - 1].pairVec, reducedVec[dim].tid, reducedVec[dim - 1].tid,
            dimFactors, splitComplex[dim - 1], sortedSplitFilt[dim - 1], sortedSplitFilt[dim], twoBelowCells, dim);
    }
}

vector<MorseReducedResult> performPersistence(const Complex& comp)
{
    /*
        Actually perform the persistence algorithm:
        First reduce all of the boundary operators to LU normal form.
        Then extract barcodes and cycle reps from the reduced complex.
    */

    const vector<SparseBool>& splitComplex          = comp.splitCells;
    const vector<uint32_t>& dimPatt                 = comp.dimPatt;
    const vector<vector<uint32_t>>& sortedSplitFilt = comp.splitFiltInds;

    const uint32_t numDim = numberOfDim(dimPatt);

    vector<MorseReducedResult> resultVec(numDim + 1);
    vector<vector<uint32_t>> prePairs(numDim + 1, vector<uint32_t>{});

    // spla::testSilents();

    // Calculating cohomology: start from low dimensions
    // Each iteration factorizes boundary from dim to (dim - 1)
    for (uint32_t dim = 1; dim < numDim; ++dim)
    {
        // Get the high cell indices from the pairs calculated in the lower dimension
        vector<uint32_t> lowbasisnames(resultVec[dim - 1].pairVec.size());
        std::transform(resultVec[dim - 1].pairVec.begin(), resultVec[dim - 1].pairVec.end(), lowbasisnames.begin(),
                       [](const eirene::MorsePair& pair) { return pair.highInd; });

        const SparseBool& dimBoundary = splitComplex[dim];
        uint32_t numHighCells = dimBoundary.numCols(), numLowCells = splitComplex[dim - 1].numCols();

        // sanity check
        if (numLowCells != dimBoundary.numRows) std::cerr << "Error! Wrong row/col dims!" << std::endl;

        // low lab will hold the indices of cells of dimension (dim - 1) which do not appear in
        // lowbasisnames, i.e. were unpaired in any previous iterations
        vector<uint32_t> lowlab = uniqueIntsNotInUnsortedUpTo(lowbasisnames, numLowCells);

        // also, low lab will be in permuted order
        // sortedSplitFilt is indices into a vector of unique filtration values in REVERSE
        // sorted order, so a bigger index means the cell is born earlier
        vector<float> subFilt = vector<float>(lowlab.size());
        std::transform(lowlab.begin(), lowlab.end(), subFilt.begin(),
                       [&](uint32_t ind) { return sortedSplitFilt[dim - 1][ind]; });

        vector<uint32_t> lowCellFiltPerm = getSortPerm(subFilt);

        applyPermutation(lowlab, lowCellFiltPerm);

        // Nothing to reduce. For example, when dim == numDim - 1 and the boundary operator is zero.
        if (dimBoundary.numNonZero() == 0)
        {
            resultVec[dim].tid = lowlab;

            // the reduced coboundary is zero, but make sure it's the correct size...
            resultVec[dim].LUCoBdry.colPtrs.resize(lowlab.size() + 1);
            resultVec[dim].LUCoBdry.numRows = dimBoundary.numCols();

            continue;
        }

        vector<uint32_t> higlab(numHighCells);
        std::iota(higlab.begin(), higlab.end(), 0);

        // Transpose from boundary to coboundary. Note this also reindexes the
        // input's rows/output's columns according to lowlab
        auto dimCoBdry = spla::transposeSparseSubmatrix(dimBoundary, lowlab, higlab);

        // Bunch of copies. Good spot for some memory mgmnt..
        auto highFiltTemp = vector<uint32_t>(higlab.size());
        for (uint32_t i = 0; i < higlab.size(); ++i) highFiltTemp[i] = sortedSplitFilt[dim][higlab[i]];

        auto lowFiltTemp = vector<uint32_t>(lowlab.size());
        for (uint32_t i = 0; i < lowlab.size(); ++i) lowFiltTemp[i] = sortedSplitFilt[dim - 1][lowlab[i]];

        PairVec pairs;

        /*
            From orig eirene:
                - 	the output array tid (sic. reindexed tlab) has size equal to the number of columns of the
                    input array, NOT the column-rank of the input array
                -	the columns of M must be ordered by grain (ascending) (simplify implementation of this?)
                -	the first (rank of M) elements of tlab index the complete set
                    of nonzero columns in the reduced matrix
        */
        auto [schurCoBdry, schurPairs, tlab] = spla::morseLU(dimCoBdry, pairs, highFiltTemp, lowFiltTemp);

        std::transform(schurPairs.begin(), schurPairs.end(), schurPairs.begin(),
                       [&lowlab, &higlab](const eirene::MorsePair& pair) {
                           return eirene::MorsePair{lowlab[pair.lowInd], higlab[pair.highInd]};
                       });

        schurCoBdry.rowVals = reindex(schurCoBdry.rowVals, lowlab);

        resultVec[dim] = eirene::MorseReducedResult{schurCoBdry,             // trv and tcp
                                                    reindex(tlab, lowlab),   // tid
                                                    schurPairs,              // plo/hi
                                                    vector<vector<uint32_t>>(0)};
    }

    extractHomology(resultVec, splitComplex, dimPatt, sortedSplitFilt);

    return resultVec;
}

}   // namespace eirene
