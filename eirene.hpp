#include <variant>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <tuple>


template <typename T>
inline void print_vector(const std::vector<T> vec)
{
    if (vec.size() == 0)
    {
        std::cerr << "[]\n";
        return;
    }
    std::cerr << "[";

    for (uint32_t i = 0; i < vec.size() - 1; ++i) std::cerr << vec[i] << ", ";

    std::cerr << vec.back() << "]\n";
}


// Max: 2^32 cells
namespace eirene
{
struct EireneResult
{
    bool result = true;
};

inline std::vector<uint32_t> computeDimensionPattern(const std::vector<uint32_t>& dimVec)
{
    /*
        Transforms the dimension vector into the dimension pattern:
        Dimension vector has number of cells with dimension i in entry i
        Dimension pattern has 1 + (# of cells with dimension < i) in entry i
    */
    auto dimPatt = std::vector<uint32_t>(dimVec.size() + 1, 0);
    dimPatt[0]   = 1;
    for (uint32_t i = 0; i < dimVec.size(); ++i) dimPatt[i + 1] = dimPatt[i] + dimVec[i];
    return dimPatt;
}

struct SparseBool
{
    std::vector<uint32_t> rowVals;
    std::vector<uint32_t> colPtrs;
    size_t numRows;

    // Co or contra variance
    const enum var_t { none, co, contra } type;

    SparseBool() : SparseBool(0, 0) {}

    SparseBool(size_t numRows, size_t numCols, var_t var = var_t::none) : numRows{numRows}, type{var}
    {
        colPtrs.resize(numCols + 1);
    }

    SparseBool(std::vector<uint32_t> rowVals, std::vector<uint32_t> colPtrs, size_t numRows, var_t var = var_t::none)
        : rowVals{rowVals}, colPtrs{colPtrs}, numRows{numRows}, type{var}
    {
    }

    SparseBool(const SparseBool& other)
        : rowVals{other.rowVals}, colPtrs{other.colPtrs}, numRows{other.numRows}, type{other.type}
    {
    }

    SparseBool& operator=(SparseBool other)
    {
        std::swap(rowVals, other.rowVals);
        std::swap(colPtrs, other.colPtrs);
        std::swap(numRows, other.numRows);
        return *this;
    }

    uint32_t highDimNum_(var_t varType) const
    {
#ifndef NOERR
        if (type == var_t::none)
        {
            std::cerr << "Didn't assign variance, cant call dimnum\n";
            std::exit(2);
        }
#endif
        // defaults to contravariant if it's none
        return varType == var_t::co ? numCols() : numRows;
    }

    uint32_t highDimNum() const { return highDimNum_(type); }

    uint32_t lowDimNun() const
    {
        var_t newType = (type == var_t::none) ? type : (type == var_t::co ? var_t::contra : var_t::co);
        return highDimNum_(newType);
    }

    uint32_t numCols() const { return colPtrs.size() - 1; }

    uint32_t numNonZero() const { return rowVals.size(); }

    void emplace_at(const bool toPlace, const uint32_t rowInd, const uint32_t colInd)
    {
        // Stitched in, not tested
        const uint32_t startInd = colPtrs[colInd], stopInd = colPtrs[colInd + 1];

        // If the column in nonempty, see if the location is already nonzero
        for (uint32_t i = startInd; i < stopInd && rowVals[i] <= rowInd; ++i)
        {
            // if it is, either erase it or overwrite it, depending
            if (rowVals[i] == rowInd)
            {
                if (!toPlace)
                {
                    rowVals.erase(rowVals.begin() + i);

                    for (i = colInd + 1; i < colPtrs.size(); ++i) this->colPtrs[i] -= 1;

                    return;
                }

                return;
            }
        }

        // if we're here, the location is zero
        if (!toPlace) return;

        // keep rowInds in increasing order for each column
        uint32_t insertInd = startInd;
        while (insertInd < stopInd && rowVals[insertInd] < rowInd) insertInd++;

        rowVals.insert(rowVals.begin() + insertInd, rowInd);

        for (uint32_t i = colInd + 1; i < colPtrs.size(); ++i) colPtrs[i]++;
    }
};

inline std::ostream& operator<<(std::ostream& out, const SparseBool& obj)
{
    const uint32_t numRows = obj.numRows, numCols = obj.numCols();
    uint32_t toWrite = 0;

    for (uint32_t irow = 0; irow < numRows; ++irow)
    {
        for (uint32_t icol = 0; icol < numCols; ++icol)
        {
            uint32_t colStart = obj.colPtrs[icol], colEnd = obj.colPtrs[icol + 1];
            toWrite = 0;

            while (colStart < colEnd)
            {
                if (obj.rowVals[colStart] == irow)
                {
                    toWrite = 1;
                    break;
                }

                ++colStart;
            }

            out << toWrite << " ";
        }
        out << "\n";
    }

    return out;
}

struct MorsePair
{
    uint32_t lowInd;
    uint32_t highInd;
};

inline std::ostream& operator<<(std::ostream& out, const MorsePair& pair)
{
    out << "(" << pair.lowInd << "," << pair.highInd << ") ";
    return out;
}

struct PairVec : std::vector<MorsePair>
{
   public:
    PairVec(uint32_t sz = 0) : std::vector<MorsePair>(sz) {}
    void set_marker(uint32_t ind) { marker = ind; }
    const auto get_marker() const { return this->begin() + marker; }
    void advance_marker(uint32_t del) { marker += del; }

   private:
    uint32_t marker;
};

struct MorseReducedResult
{
    SparseBool LUCoBdry;   // trv and tcp
    std::vector<uint32_t> tid;
    PairVec pairVec;   // phi and plo
    std::vector<std::vector<uint32_t>> cycleReps;
};

struct Complex
{
    // A persistence-filtered (cellular) chain complex
    std::vector<SparseBool> splitCells;
    std::vector<uint32_t> dimPatt;

    std::vector<std::vector<uint32_t>> splitFiltInds;
    std::vector<float> uniqFiltVals;

    // Lazy<MorseReducedResult> mrr;
};

using MaybeReduced = std::variant<std::vector<MorseReducedResult>, std::string>;

struct SubComplex
{
    // A subset of an abstract complex is closed in the Alexandrov
    // topology iff it's a subcomplex i.e. downward closed
    SubComplex(const Complex& comp, std::vector<uint32_t> cells, const std::vector<uint32_t>& eulerVec)
        : cells{std::move(cells)}, dimPatt{computeDimensionPattern(eulerVec)}
    {
        (void)comp;
    }

    // Must be sorted by increasing cell dim/index
    std::vector<uint32_t> cells, dimPatt;

    bool is_closed() const;
};

using IndVec = std::vector<uint32_t>;

std::variant<eirene::Complex, std::string> toComplex(const SparseBool& adj, const std::vector<float>& filts,
                                                     const std::vector<uint32_t>& cellDims, const bool = true);

/*
struct SubComplex
{
     *  Represents a quotient complex of the whole cell complex.
     *  Specifically, the whole relative to the subset.
    std::vector<MorseReducedResult> normalForm;
    std::vector<uint32_t> subComplexCells;

   private:
    Complex whole;
};
*/

MaybeReduced eirene(const SparseBool& chainComplex, const std::vector<float>& filtVals,
                    const std::vector<uint32_t>& dimVec);

}   // namespace eirene
