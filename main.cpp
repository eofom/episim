#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <array>
#include <memory>

#define OUTSTR(x) (std::string(#x": ") + std::to_string(x))

class RandomGenerator {
public:
    // singleton pattern
    RandomGenerator(RandomGenerator const&) = delete;
    RandomGenerator(RandomGenerator &&) = delete;
    RandomGenerator operator=(RandomGenerator const&) = delete;
    RandomGenerator operator=(RandomGenerator &&) = delete;

    static RandomGenerator& getInstance(int seed = 0) {
        static RandomGenerator instance(seed);
        return instance;
    }

    static void setSeed(int seed) {
        auto &instance = getInstance();
        instance.generator = std::mt19937(seed);
    }

    // generate random int in [from, to)
    int uniformInt(int from, int to) {
        // -1 is here because uniform_int_distribution generates from a to b inclusevly
        return std::uniform_int_distribution<int>(from, to - 1)(generator);
    }

    int uniformIntWithExclusion(int from, int to, int exclusion) {
        int omega = uniformInt(from, to - 1);
        return omega == exclusion ? to - 1 : omega;
    }

    int uniformIntUpTo(int to) {
        return uniformInt(0, to);
    }

    int uniformIntUpToWithExclusion(int to, int exclusion) {
        return uniformIntWithExclusion(0, to, exclusion);
    }

    int randomChoice(const std::vector<double> &probabilities) {
        double cumulativeProbability = 0;
        double omega = realDistribution(generator);
        for (int i = 0; i < probabilities.size(); ++i) {
            cumulativeProbability += probabilities.at(i);
            if (omega < cumulativeProbability) {
                return i;
            }
        }
        return probabilities.size();
    }

    bool randomEvent(double probability) {
        return randomChoice({probability}) == 0;
    }

private:
    explicit RandomGenerator(int seed) : realDistribution(0., 1.0), generator(seed) {}

    std::uniform_real_distribution<double> realDistribution;
    std::mt19937 generator;
};

RandomGenerator &getRandomGenerator() {
    return RandomGenerator::getInstance();
}

void seedRandomGenerator(int seed) {
    RandomGenerator::setSeed(seed);
}

struct Disease {
    Disease(double avgDuration, double deathRate, double transmissionProbablity = 0.)
    : transmissionProbability(transmissionProbablity) {
        dailyPersistenceProbability = double(avgDuration - 1.) / avgDuration;
        dailyCureProbability = (1. - dailyPersistenceProbability) * (1. - deathRate);
    }

    double dailyPersistenceProbability;
    double dailyCureProbability;
    double transmissionProbability;

    [[nodiscard]] double deathDailyProbability() const {
        return 1 - dailyPersistenceProbability - dailyCureProbability;
    }
};

enum class State {
    kHealthy = 0,
    kIll,
    kImmune,
    kDead
};

class Person {
public:
    Person(State state) : state_(state), disease_(nullptr) {}
    Person(std::shared_ptr<Disease> disease) : state_(State::kIll), disease_(std::move(disease)) {}

    void iterateState() {
        if (disease_ == nullptr) {
            return;
        }
        if (state_ == State::kHealthy) {
            state_ = State::kIll;
            return;
        }
        assert(state_ == State::kIll);
        auto &randomGenerator = getRandomGenerator();
        auto outcome = randomGenerator.randomChoice({disease_->dailyPersistenceProbability, disease_->dailyCureProbability});
        if (outcome == 1) {
            disease_ = nullptr;
            state_ = State::kImmune;
        } else if (outcome == 2) {
            disease_ = nullptr;
            state_ = State::kDead;
        }
    }

    void meet(Person &other) {
        tryTransmitTo(other);
        other.tryTransmitTo(*this);
    }

    [[nodiscard]] State getState() const {
        return state_;
    }

    bool noDisease() {
        return disease_ == nullptr;
    }

private:
    void tryTransmitTo(Person &other) {
        if (state_ != State::kIll || other.state_ != State::kHealthy) {
            return;
        }
        assert(disease_ != nullptr);
        if (getRandomGenerator().randomEvent(disease_->transmissionProbability)) {
            other.disease_ = disease_;
        }
    }

    State state_;
    std::shared_ptr<Disease> disease_;
};

struct IllnessResult {
    int duration;
    State state;
};

IllnessResult simulateIllness(const Disease &disease) {
    Person person(std::make_shared<Disease>(disease));
    std::weak_ptr<Disease> diseaseWatcher;
    int duration{};
    while (person.getState() == State::kIll) {
        ++duration;
        person.iterateState();
    }
    return {duration, person.getState()};
}

struct IllnessBatchResult {
    double avgDuration;
    double deathRate;
};

IllnessBatchResult simulateIllnessBatch(const Disease &disease, int iterations) {
    int64_t totalDuration = 0;
    int deathCount = 0;

    for (int i = 0; i < iterations; ++i) {
        auto results = simulateIllness(disease);
        totalDuration += results.duration;
        deathCount += (results.state == State::kDead);
    }
    return {totalDuration / static_cast<double>(iterations), deathCount / static_cast<double>(iterations)};
}

struct EpidemicResult {
    int deadCount;
    int immuneCount;
    int healthyCount;
    int duration;

    void print() {
        std::cout << OUTSTR(deadCount) << " " << OUTSTR(immuneCount) << " "
                  << OUTSTR(healthyCount) << " " << OUTSTR(duration) << std::endl;
    }
};

EpidemicResult simulateEpidemic(const Disease &disease, int populationSize, int initiallyInfected,
                                double avgContactsPerPerson, bool verbose = false) {
    assert(populationSize >= initiallyInfected);
    assert(initiallyInfected >= 0);
    assert(avgContactsPerPerson <= populationSize - 1);
    if (initiallyInfected == 0) {
        return {0, 0, populationSize, 0};
    }

    auto patientZero = std::make_shared<Disease>(disease); // used only for initialization
    std::vector<Person> population(initiallyInfected, patientZero);
    population.resize(populationSize, State::kHealthy);
    assert(population.back().noDisease());
    std::weak_ptr<Disease> diseaseWatcher = patientZero;
    patientZero = nullptr;

    auto meetings_per_day = static_cast<int>(populationSize * avgContactsPerPerson);

    int duration = 0;
    while (!diseaseWatcher.expired()) {
        ++duration;
        for (int i = 0; i < meetings_per_day; ++i) {
            int firstIndex = getRandomGenerator().uniformIntUpTo(population.size());
            int secondIndex = getRandomGenerator().uniformIntUpToWithExclusion(population.size(), firstIndex);
            population.at(firstIndex).meet(population.at(secondIndex));
        }

        int illCount{};
        for (auto &person : population) {
            person.iterateState();
            illCount += person.getState() == State::kIll;
        }
        if (verbose) {
            std::cout << "day " << duration << " ill: " << illCount << std::endl;
        }
    }

    EpidemicResult results{};
    results.duration = duration;
    for (const auto &person : population) {
        results.deadCount += person.getState() == State::kDead;
        results.healthyCount += person.getState() == State::kHealthy;
        results.immuneCount += person.getState() == State::kImmune;
        assert(person.getState() != State::kIll);
    }
    return results;
}

bool near(double lhs, double rhs, double margin) {
    return std::abs(lhs - rhs) < margin;
}

void testDiseaseModel() {
    double avgDuration = 10.;
    double deathRate = 0.1;
    Disease disease(avgDuration, deathRate);
    auto results = simulateIllnessBatch(disease, 1e5);
    assert(near(avgDuration, results.avgDuration, 0.1));
    assert(near(deathRate, results.deathRate, 0.1));
}

void testRandomEvents(double probability, double margin) {
    double iterations = 1000000;
    double events = 0;
    for (int i = 0; i < iterations; ++i) {
        if (getRandomGenerator().randomEvent(probability)) {
            ++events;
        }
    }
    assert(near(events / iterations, probability, margin));
}

void testZeroProbabilityEvent() {
    testRandomEvents(0., 1e-15);
}

void testEvent() {
    testRandomEvents(0.3, 1e-3);
}

template<typename Test>
void runTest(Test test, const std::string &description) {
    static int test_number = 1;
    std::cout << "Test " << test_number++ << ": " << description;
    test();
    std::cout << " OK" << std::endl;
}

void runAllTests() {
    seedRandomGenerator(0);
    runTest(testDiseaseModel, "disease is correctly constructed via avg duration and death rate");
    runTest(testZeroProbabilityEvent, "zero probability event never occurs");
    runTest(testEvent, "test event has proper likelihood");
}

int main() {
//    runAllTests();
    seedRandomGenerator(std::chrono::system_clock::now().time_since_epoch().count());
    Disease disease(10, 0.02, 0.05);
    auto results = simulateEpidemic(disease, 1e5, 100, 2., true);
    results.print();

    return 0;
}
