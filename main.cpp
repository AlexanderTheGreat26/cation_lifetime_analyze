#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <tuple>
#include <sstream>
#include <algorithm>
#include <utility>
#include <array>
#include <memory>
#include <stdexcept>
#include <cmath>


typedef std::vector<std::string> frame;
typedef std::vector<frame> frames;


const int step = 100;
const double dt = 0.5;


frames read (const std::string & name, const std::vector<std::string> & cation_hypostasis);

void data_file_creation (const std::string & name, frames & data, const int & step);

std::vector<int> uniq_lifes (frames & cation_frames);

std::vector<std::pair<double, int>>
life_histogram_creation (std::vector<int> &cation_times, const double & dt, const double & step);

void file_creation (const std::string & file_name, std::vector<std::pair<double, int>> & data);


int main() {

    auto Zundels_lines = read("CoMs.pos", {"H5O2\t"});
    //data_file_creation("Zundels.only", Zundels_lines, step);

    auto H3O_lines = read("CoMs.pos", {"H3O\t", "H6O2\t", "H9O4\t"});
    //data_file_creation("H3O.only", H3O_lines, step);

    std::vector<int> Zundel_times = std::move(uniq_lifes(Zundels_lines));
    std::vector<int> H3O_times = std::move(uniq_lifes(H3O_lines));
    // And the last step - histogram via GNUPlot.

    auto uniq_Z = life_histogram_creation(Zundel_times, dt, step);
    file_creation("Zundels_uniq", uniq_Z);

    auto uniq_H3O = life_histogram_creation(H3O_times, dt, step);
    file_creation("H3O_uniq", uniq_H3O);

    return 0;
}


bool contain (const std::string & word, const std::string & line) {
    return line.find(word) != std::string::npos;
}


frames read (const std::string & name, const std::vector<std::string> & cation_hypostasis) {
    frames cations;
    std::string line;
    int i = 0;
    std::ifstream fin(name);
    if (!fin.is_open()) throw std::runtime_error("Error opening file.");
    while (!fin.eof())
        while (getline(fin, line)) {
            bool timestep_line = contain("Timestep", line);
            if (!(timestep_line || contain("ID", line))) {

                for (const auto & cation: cation_hypostasis)
                    if (contain(cation, line) && !contain("Cl", line))
                        cations[i].emplace_back(line);

            } else if (timestep_line && !cations.empty()) {
                ++i;
            } else cations.resize(cations.size() + 1);
        }
    return cations;
}


void data_file_creation (const std::string & name, frames & data, const int & step) {
    std::ofstream fout;
    fout.open(name, std::ios::trunc);
    for (int i = 0; i < data.size()-1; ++i) {
        fout << "#\n# Timestep\t" << i*step << '\n';
        for (const auto & line: data[i])
            fout << line << '\n';
    }
    fout.close();
}


template<typename T>
T fromString(const std::string& s){
    std::istringstream iss(s);
    T res;
    iss >> res;
    return res;
}


bool any_of (std::string & current_cation, frame & cation_frame, int & index) {
    for (int i = 0; i < cation_frame.size(); ++i)
        if (fromString<int>(current_cation) == fromString<int>(cation_frame[i])) {
            index = i;
            return true;
        }
    return false;
}


bool is_dead (std::string & uniq_cation, int & frame_number, frames & cation_frames, const int & life_step) {
    int cation_in_frame_pos;
    for (int i = frame_number; i < frame_number + life_step; ++i)
        if (any_of(uniq_cation, cation_frames[i], cation_in_frame_pos)) {
            cation_frames[i].erase(cation_frames[i].begin() + cation_in_frame_pos);
            cation_frames[frame_number].emplace_back(uniq_cation);
            return false;
        }
    return true;
}


int alive_for (std::string & uniq_cation, frames & cation_frames, int & frame_number) {
    int cation_in_frame_pos, lifetime = 0;
    int current_step, max_life_step = 4; // Empirical parameter.
    for (int i = frame_number; i < cation_frames.size(); ++i) {
        if (any_of(uniq_cation, cation_frames[i], cation_in_frame_pos)) {
            cation_frames[i].erase(cation_frames[i].begin()+cation_in_frame_pos);
            ++lifetime;
            current_step = i;
            continue;
        } else if (i + max_life_step > cation_frames.size()) {
            max_life_step = cation_frames.size() - i;
        } if (i <= current_step + max_life_step &&
              !is_dead(uniq_cation, i, cation_frames, max_life_step)) {
            continue;
        } else break;
    }
    return lifetime;
}


std::vector<int> uniq_lifes (frames & cation_frames) {
    std::vector<int> result;
    std::string cation;
    for (int i = 0; i < cation_frames.size(); ++i) // 2000 timesteps for relaxation.
        while (!cation_frames[i].empty())
            for (int j = 0; j < cation_frames[i].size(); ++j)
                result.emplace_back(alive_for(cation_frames[i][j], cation_frames, i));
    return result;
}


std::vector<std::pair<double, int>>
groups (std::vector<int> & borders, std::vector<int> & data, const double & dt, const int & step) {
    std::vector<std::pair<double, int>> result (borders.size()-1);
    for (int & j : data)
        for (int i = 1; i < borders.size()-1; ++i)
            if (j >= borders[i-1] && j < borders[i])
                ++result[i-1].second;
    for (int i = 1; i < result.size(); ++i)
        result[i-1].first = double(i) * double(step) * dt;
    return result;
}



std::vector<std::pair<double, int>>
life_histogram_creation (std::vector<int> &cation_times, const double & dt, const double & step) {
    //int first_group_border = 1;
    int shortest_life = *std::min_element(cation_times.begin(), cation_times.end());
    int longest_life = *std::max_element(cation_times.begin(), cation_times.end());
    std::vector<int> group_borders(longest_life-shortest_life+1);
    std::generate(group_borders.begin(), group_borders.end(), [&] {return shortest_life++;});
    std::vector<std::pair<double, int>> histogram = groups(group_borders, cation_times, dt, step);
    //file_creation(cation_name, histogram);
    return histogram;
}


void file_creation (const std::string & file_name, std::vector<std::pair<double, int>> & data) {
    std::ofstream fout;
    fout.open(file_name, std::ios::trunc);
    for (auto & i : data)
        fout << i.first << '\t' << i.second << '\n';
    fout.close();
}