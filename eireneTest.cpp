#include "eirene.hpp"

using std::vector;

auto testLine()
{
    eirene::SparseBool testComplex(3, 3);
    testComplex.emplace_at(true, 0, 2);
    testComplex.emplace_at(true, 1, 2);

    std::cout << testComplex << "\n\n";

    std::vector<float> testFilt{0.01, 0.02, 0.03};
    std::vector<uint32_t> eulerVec{2, 1};

    return eirene::eirene(testComplex, testFilt, eulerVec);
}

auto testTriangle()
{
    eirene::SparseBool testComplex(6, 6);

    testComplex.emplace_at(true, 0, 3);
    testComplex.emplace_at(true, 1, 3);

    testComplex.emplace_at(true, 1, 4);
    testComplex.emplace_at(true, 2, 4);

    testComplex.emplace_at(true, 0, 5);
    testComplex.emplace_at(true, 2, 5);

    std::vector<float> testFilt{0.01, 0.02, 0.03, 0.04, 0.05, 0.06};
    std::vector<uint32_t> eulerVec{3, 3};

    //    std::cerr << testComplex << std::endl;

    return eirene::eirene(testComplex, testFilt, eulerVec);
}

auto testSphere()
{
    eirene::SparseBool testComplex(6, 6);

    testComplex.emplace_at(true, 0, 2);
    testComplex.emplace_at(true, 1, 2);

    testComplex.emplace_at(true, 0, 3);
    testComplex.emplace_at(true, 1, 3);

    testComplex.emplace_at(true, 2, 4);
    testComplex.emplace_at(true, 3, 4);

    testComplex.emplace_at(true, 2, 5);
    testComplex.emplace_at(true, 3, 5);

    std::cerr << testComplex << std::endl;

    auto testFilt = std::vector<float>{0.01, 0.02, 0.03, 0.04, 0.05, 0.06};
    auto eulerVec = std::vector<uint32_t>{2, 2, 2};

    return eirene::eirene(testComplex, testFilt, eulerVec);
}

auto testLingareddy()
{
    /*
     * Random complex, no filtration, from https://math.uchicago.edu/~may/REU2018/REUPapers/Lingareddy.pdf.
     */
    eirene::SparseBool testComplex(13, 13);

    testComplex.emplace_at(true, 0, 5);
    testComplex.emplace_at(true, 0, 6);
    testComplex.emplace_at(true, 0, 7);

    testComplex.emplace_at(true, 1, 5);
    testComplex.emplace_at(true, 1, 8);

    testComplex.emplace_at(true, 2, 6);
    testComplex.emplace_at(true, 2, 9);

    testComplex.emplace_at(true, 3, 10);
    testComplex.emplace_at(true, 3, 9);
    testComplex.emplace_at(true, 3, 8);
    testComplex.emplace_at(true, 3, 7);
    
    testComplex.emplace_at(true, 4, 10);

    testComplex.emplace_at(true, 5, 11);
    testComplex.emplace_at(true, 7, 11);
    testComplex.emplace_at(true, 8, 11);

    testComplex.emplace_at(true, 6, 12);
    testComplex.emplace_at(true, 7, 12);
    testComplex.emplace_at(true, 9, 12);

    std::cerr << testComplex << std::endl;

    auto testFilt = std::vector<float>(13, 0.01);
    auto eulerVec = std::vector<uint32_t>{5, 6, 2};

    return eirene::eirene(testComplex, testFilt, eulerVec);
}

int main()
{
    // testLingareddy();
    auto out = std::get<std::vector<eirene::MorseReducedResult>>(testLine());
    uint32_t dim = 0;

    for (auto& homInfo : out)
    {
        std::cerr << "Cycle reps from dim " << dim << "\n";

        for (auto& repInfo : homInfo.cycleReps)
            print_vector<>(repInfo);

        ++dim;
    }
    return 0;
}
