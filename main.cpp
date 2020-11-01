#include <iostream>
#include <random>
#include <chrono>


class Simulator {
public:
    enum class State {
        kHealthy = 0,
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

    Simulator(int average_illness_length, double death_rate, int random_precision)
    : generator(std::chrono::system_clock::now().time_since_epoch().count()), distribution(1, random_precision) {
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
            state = ProgressIllness();
        }
        return {duration, state};
    }

    SimulationResult SimulateBatch(int count) {
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

        return {static_cast<double>(total_duration) / (count - death_count), static_cast<double>(death_count) / count};
    }

private:
    State ProgressIllness() {
        int random_token = distribution(generator);
        if (random_token <= ill_cutoff) {
            return State::kIll;
        } else if (random_token < get_well_cutoff) {
            return State::kImmune;
        } else {
            return State::kDead;
        }
    }

    std::mt19937 generator;
    std::uniform_int_distribution<int> distribution;
    int ill_cutoff;
    int get_well_cutoff;
};

int main() {
    auto begin = std::chrono::steady_clock::now();
    srand(std::chrono::system_clock::now().time_since_epoch().count());
    Simulator simulator(10, 0.3, 1e6);
    auto simulation_result = simulator.SimulateBatch(1e7);
    std::cout << "avg length: " << simulation_result.avg_ill_time << " death rate: " << simulation_result.death_rate << std::endl;
    auto end = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms" << std::endl;
    return 0;
}
