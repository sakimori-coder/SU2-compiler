#pragma once

#include <chrono>
#include <unordered_map>
#include <string>
#include <iostream>
#include <iomanip>
#include <mutex>

namespace su2compiler::util {

using clock_t = std::chrono::steady_clock;
using micro   = std::chrono::duration<double,std::micro>;

/* ================================================================== */
class Timer {
    clock_t::time_point t0_{};
    micro               acc_{0};
    bool                running_ = false;
public:
    void   start() noexcept { if(!running_){ t0_=clock_t::now(); running_=true; } }
    void   stop()  noexcept { if(running_){ acc_ += clock_t::now()-t0_; running_=false; } }
    void   reset() noexcept { acc_ = micro{0}; running_ = false; }
    micro  elapsed() const noexcept {
        if(running_)
            return acc_ + std::chrono::duration_cast<micro>(clock_t::now() - t0_);
        return acc_;
    }
};



/* ================================================================== */
class Profiler {
    std::unordered_map<std::string,Timer> table_;
    std::mutex mtx_;
    bool auto_report_ = true;                 // toggle for destructor output
    clock_t::time_point total_start = clock_t::now();
public:
    static Profiler& instance(){ static Profiler p; return p; }
    /* get named timer */
    Timer& operator[](const std::string& name){ return table_[name]; }

    // manual start/stop so callers can time across multiple scopes
    void start(const std::string& name){ table_[name].start(); }
    void stop (const std::string& name){ table_[name].stop (); }

    /* pretty print summary */
    void report(std::ostream& os = std::cerr) const {
        if(table_.empty()) return;
        os << "\n──── timing summary ─────────────────────────────\n";
        double total = std::chrono::duration_cast<micro>(clock_t::now() - total_start).count();
        for(auto& [k,tm]:table_){
            double us = tm.elapsed().count();
            os << std::setw(24) << k << " : "
               << std::setw(10) << std::fixed << std::setprecision(3) << us/1000.0
               << " ms  (" << std::setprecision(2) << (total? (us/total*100.0):0) << "%)\n";
        }
    }
    void set_auto_report(bool v){ auto_report_ = v; }

    /* automatic summary at program exit */
    ~Profiler(){ if(auto_report_) report(); }

    /* helper RAII object */
    class Section {
        Timer& t_;
    public:
        Section(const std::string& name) : t_(Profiler::instance()[name]){ t_.start(); }
        ~Section(){ t_.stop(); }
    };
};

}