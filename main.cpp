#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <array>

#define OUTSTR(x) (std::string(#x": ") + std::to_string(x))

class Simulator {
public:
    enum class State {
        kHealthy = 0,
        kInIncubationPeriod,
        kIll,
        kImmune,
        kDead
    };

    struct IllnessResult {
        int duration;
        State state;
    };

    struct SimulationResult {
        double avg_ill_time;
        double death_rate;
    };

    struct EpidemicResult {
        int dead_count;
        int immune_count;
        int healthy_count;
        int duration;
    };

    Simulator(int average_illness_length, double death_rate, int random_precision)
    : generator(std::chrono::system_clock::now().time_since_epoch().count()), random_precision(random_precision) {
        double continue_being_ill_probablility = double(average_illness_length - 1.) / average_illness_length;
        double get_well_probability = (1. - continue_being_ill_probablility) * (1. - death_rate);

        ill_cutoff = static_cast<int>(continue_being_ill_probablility * random_precision);
        get_well_cutoff = static_cast<int>((continue_being_ill_probablility + get_well_probability) * random_precision);
    }

    IllnessResult SimulateIllness() {
        int duration = 0;
        State state = State::kIll;
        while (state == State::kIll) {
            ++duration;
            state = ProgressIllness(state);
        }
        return {duration, state};
    }

    SimulationResult SimulateIllnessBatch(int count) {
        int death_count = 0;
        int64_t total_duration = 0;
        for (int i = 0; i < count; ++i) {
            auto illness_result = SimulateIllness();
            if (illness_result.state == State::kDead) {
                ++death_count;
            } else {
                total_duration += illness_result.duration;
            }
        }

        return {static_cast<double>(total_duration) / (count - death_count),
                static_cast<double>(death_count) / count};
    }

    // TODO: simulate many
    EpidemicResult SimulateEpidemic(const int population_size, const int initially_infected,
                                    const double average_contact_per_person, const double infection_probability) {
        assert(population_size >= initially_infected);
        assert(average_contact_per_person <= population_size - 1);
        std::vector<State> population(initially_infected, State::kIll);
        population.resize(population_size, State::kHealthy);

        auto meetings_per_day = static_cast<int>(static_cast<double>(population_size) * average_contact_per_person);
        int ill_count = initially_infected;
        int dead_count = 0;
        int immune_count = 0;
        int days = 0;
        while (ill_count > 0) {
            ++days;

            for (int i = 0; i < meetings_per_day; ++i) {
                auto first_person_index = std::uniform_int_distribution<int>(0, population.size() - 2)(generator);
                auto &first_person = population.at(first_person_index);
                auto second_person_index = std::uniform_int_distribution<int>(first_person_index + 1, population_size - 1)(generator);
                auto &second_person = population.at(second_person_index);

                if (first_person == State::kIll && second_person == State::kHealthy) {
                    if (RandomEvent(infection_probability)) {
                        second_person = State::kInIncubationPeriod;
                    }
                } else if (second_person == State::kIll && first_person == State::kHealthy) {
                    if (RandomEvent(infection_probability)) {
                        first_person = State::kInIncubationPeriod;
                    }
                }
            }
            for (auto &person : population) {
                auto person_before = person;
                person = ProgressIllness(person);
                if (person != person_before) {
                    if (person == State::kIll) {
                        ++ill_count;
                    } else if (person_before == State::kIll) {
                        --ill_count;
                        if (person == State::kDead) {
                            ++dead_count;
                        } else if (person == State::kImmune) {
                            ++immune_count;
                        }
                    }
                }
            }
        }

        return {dead_count, immune_count, population_size - dead_count - immune_count, days};
    }

private:
    State ProgressIllness(State state) {
        switch (state) {
            case State::kInIncubationPeriod:
                return State::kIll;
            case State::kIll:
                // TODO: check performance for static std::distribution
                if (int random_token = RandomToken(); random_token < ill_cutoff) {
                    return State::kIll;
                } else if (random_token < get_well_cutoff) {
                    return State::kImmune;
                } else {
                    return State::kDead;
                }
            default:
                return state;
        }
    }

    bool RandomEvent(const double probability) {
        int random_token = RandomToken();
        return random_token < random_precision * probability;
    }

    int RandomToken() {
        return std::uniform_int_distribution<int>(1, random_precision)(generator);
    }

    std::mt19937 generator;
    int random_precision;
    int ill_cutoff;
    int get_well_cutoff;
};

void TestDeathRateAndDuration() {
    // TODO: rewrite as a test
    // TODO: add more tests
    Simulator simulator(10, 0.3, 1e6);
    auto simulation_result = simulator.SimulateIllnessBatch(1e6);
    std::cout << OUTSTR(simulation_result.avg_ill_time) << " " << OUTSTR(simulation_result.death_rate) << std::endl;
}

int main() {
    auto begin = std::chrono::steady_clock::now();
    srand(std::chrono::system_clock::now().time_since_epoch().count());

    Simulator simulator(20, 0.03, 1e6);
    auto epidemic_result = simulator.SimulateEpidemic(10000, 100, 3, 0.01);
    std::cout << OUTSTR(epidemic_result.dead_count) << "\n" << OUTSTR(epidemic_result.immune_count) << "\n"
              << OUTSTR(epidemic_result.healthy_count) << "\n" << OUTSTR(epidemic_result.duration) << std::endl;
    auto end = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms" << std::endl;
    return 0;
}
