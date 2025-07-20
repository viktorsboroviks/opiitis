#include <boost/json/src.hpp>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include "aviize.hpp"
#include "garaaza.hpp"
#include "iestaade.hpp"

namespace opiitis {

// forward declarations
// clang-format off
struct Settings;
template <typename T> struct Runtime;
template <typename T> struct StateMachine;
// clang-format on

// all classes packaged as a single template parameter
struct Context {
    using Settings     = Settings;
    using Runtime      = Runtime<Context>;
    using StateMachine = StateMachine<Context>;
};

class Schedule {
private:
    struct Period {
        size_t i_start;
        size_t period;
    };

    std::vector<Period> _periods;

    void _init_periods(const std::string &config_path,
                       const std::string &key_path,
                       std::vector<Period> &periods)
    {
        const boost::json::array periods_json =
                iestaade::value_from_json(config_path, key_path).as_array();
        assert(periods_json.size() > 0);
        assert(periods.empty());

        for (const boost::json::value &period_json : periods_json) {
            const size_t i_start =
                    period_json.as_array()[0].to_number<size_t>();
            const size_t period =
                    period_json.as_array()[1].to_number<size_t>();
            periods.emplace_back(Period{i_start, period});
        }
    }

    static size_t _i_period(size_t i, const std::vector<Period> &periods)
    {
        assert(!periods.empty());

        // check i
        if (i < periods.front().i_start) {
            return 0;
        }

        // go over all periods
        for (size_t i_period = 0; i_period + 1 < periods.size(); i_period++) {
            // find first that fits
            if (i >= periods[i_period].i_start &&
                i < periods[i_period + 1].i_start) {
                return i_period;
            }
        }

        // only last period remains
        if (i >= periods.back().i_start) {
            return periods.size() - 1;
        }
        else {
            assert(false);
            return 0;
        }
    }

public:
    Schedule() {}

    Schedule(const std::string &config_path, const std::string &key_path)
    {
        _init_periods(config_path, key_path, _periods);
    }

    bool is_time(size_t i) const
    {
        assert(!_periods.empty());
        if (i < _periods.front().i_start) {
            return false;
        }
        const size_t i_period = _i_period(i, _periods);
        return i % _periods[i_period].period == 0;
    }
};

struct Settings {
    // sim
    size_t n_cycles = 1000;

    // progress
    size_t progress_update_period = 1;
    std::string log_file_name     = "log.csv";
    std::string stats_file_name   = "stats.txt";

    Schedule rec_schedule;

    Settings() {}

    void from_config(const std::string &config_path,
                     const std::string &key_prefix)
    {
        // clang-format off
        n_cycles                = iestaade::from_json<size_t>     (config_path, key_prefix + "n_cycles");

        progress_update_period  = iestaade::from_json<size_t>     (config_path, key_prefix + "progress_update_period");
        log_file_name           = iestaade::from_json<std::string>(config_path, key_prefix + "log_file_name");
        stats_file_name         = iestaade::from_json<std::string>(config_path, key_prefix + "stats_file_name");

        rec_schedule            = Schedule(config_path, key_prefix + "rec_periods");
        // clang-format on
    }

    Settings(const std::string &config_path, const std::string &key_prefix)
    {
        from_config(config_path, key_prefix);
    }
};

template <typename T = Context>
struct Runtime {
    // status
    size_t smfg_i  = 0;
    bool smfg_done = false;
    size_t cycle_i = 0;

    // performance
    std::chrono::time_point<std::chrono::steady_clock> start_time;
    double cycle_time_us = 0;

    // recording
    std::queue<size_t> rec_cycles_queue;
    bool do_rec = false;

    // data display and logging
    aviize::Progress progress;
    std::ofstream log_f;

    std::string stats() const
    {
        const std::chrono::time_point<std::chrono::steady_clock> time_now =
                std::chrono::steady_clock::now();
        const double runtime_s =
                std::chrono::duration_cast<std::chrono::seconds>(time_now -
                                                                 start_time)
                        .count();
        const double cycle_s         = cycle_i / runtime_s;
        const size_t first_col_width = 16;
        std::stringstream ss{};
        ss << std::left << std::setw(first_col_width) << "cycles" << cycle_i
           << std::endl;
        ss << std::left << std::setw(first_col_width) << "runtime"
           << aviize::seconds_to_hhmmss_string(runtime_s) << std::endl;
        ss << std::left << std::setw(first_col_width) << "cycle/s" << cycle_s
           << std::endl;
        return ss.str();
    }

    static std::string log_header()
    {
        std::stringstream ss;
        ss << "cycle_i";
        return ss.str();
    }

    std::string log_line() const
    {
        std::stringstream ss;
        ss << cycle_i;
        return ss.str();
    }
};

template <typename T = Context>
struct StateMachine {
    // settings
    const typename T::Settings settings;

    // runtime
    typename T::Runtime runtime;

    // functions
    typedef std::function<void(StateMachine<T> &)> sm_function_t;
    struct SMFunctionGroup {
        std::string name;  // cppcheck-suppress unusedStructMember
        bool repeat_until_done = false;
        std::vector<state_function_t> functions;
    }

    // this vector must be filled before run()
    // e.g.:
    // smfgs.push_back(SMFunctionGroup{"init", false, {init_rec_periods,
    //                                                 init_log}});
    // smfgs.push_back(SMFunctionGroup{"run", true, {decide_rec,
    //                                               select_pool,
    //                                               sample_params,
    //                                               update_state,
    //                                               rm_excess_samples}});
    std::vector<SMFunctionGroup>
            smfgs;

    explicit StateMachine(const std::string &config_path,
                          const std::string &key_prefix) :
        settings(config_path, key_prefix)
    {
    }

    void run()
    {
        runtime.start_time = std::chrono::steady_clock::now();

        assert(!smfgs.empty());

        // go over all state machine function groups (smfgs)
        // - each group can be run multiple times until smfg_done
        // - track smfg_i to know which group is running
        // - track cycle_time_us to know how long it took to run one cycle
        //   of the smfg
        for (runtime.smfg_i = 0; runtime.smfg_i < smfgs.size();
             runtime.smfg_i++) {
            const auto &smfg      = smfgs[runtime.smfg_i];
            auto cycle_begin_time = std::chrono::steady_clock::now();

            runtime.smfg_done = false;
            do {
                // go over all functions in the group
                // - abort on smfg_done
                for (const state_function_t &f : smfg.functions) {
                    if (runtime.smfg_done) {
                        break;
                    }
                    f(*this);
                }

                // update the cycle time
                const auto cycle_end_time = std::chrono::steady_clock::now();
                runtime.cycle_time_us =
                        std::chrono::duration_cast<std::chrono::microseconds>(
                                cycle_end_time - cycle_begin_time)
                                .count();
                cycle_begin_time = cycle_end_time;

                // if repeat_until_done is not set, abort
                if (!smfg.repeat_until_done) {
                    break;
                }
            } while (!runtime.smfg_done);
        }
    }
};

template <typename T>
void init_log(typename T::StateMachine &sm)
{
    auto &log_f = sm.runtime.log_f;
    assert(!log_f.is_open());
    if (sm.settings.log_file_name.empty()) {
        return;
    }

    log_f.open(sm.settings.log_file_name);
    log_f << sm.runtime.log_header();
    log_f << std::endl;
}

template <typename T>
void update_log(typename T::StateMachine &sm)
{
    auto &log_f = sm.runtime.log_f;
    assert(!sm.settings.log_file_name.empty());
    if (!log_f.is_open()) {
        return;
    }

    log_f << sm.runtime.log_line();
    log_f << std::endl;
}

template <typename T>
void progress_init(typename T::StateMachine &sm)
{
    sm.runtime.progress.n_min         = 1;
    sm.runtime.progress.n_max         = sm.settings.n_cycles;
    sm.runtime.progress.update_period = sm.settings.progress_update_period;
}

template <typename T>
void progress_text_reset(typename T::StateMachine &sm)
{
    sm.runtime.progress.reset();
}

template <typename T>
void progress_text_add_nl(typename T::StateMachine &sm)
{
    sm.runtime.progress.text += "\n";
}

template <typename T>
void progress_text_add_speed(typename T::StateMachine &sm)
{
    const double cycle_s = 1 / sm.runtime.cycle_time_us * 1000000;

    std::stringstream ss;
    ss << " n/s " << cycle_s;

    sm.runtime.progress.text += std::string(ss.str());
}

template <typename T>
void progress_text_add_total(typename T::StateMachine &sm)
{
    if (!sm.runtime.progress.last_line_empty()) {
        sm.runtime.progress.text += " ";
    }
    sm.runtime.progress.text += sm.runtime.progress.str_total(sm.cycle_i);
}

template <typename T>
void progress_text_add_pct(typename T::StateMachine &sm)
{
    if (!sm.runtime.progress.last_line_empty()) {
        sm.runtime.progress.text += " ";
    }
    sm.runtime.progress.text += sm.runtime.progress.str_pct(sm.cycle_i);
}

template <typename T>
void progress_text_add_eta(typename T::StateMachine &sm)
{
    if (!sm.runtime.progress.last_line_empty()) {
        sm.runtime.progress.text += " ";
    }
    sm.runtime.progress.text += sm.runtime.progress.str_eta(sm.cycle_i);
}

template <typename T>
void progress_print(typename T::StateMachine &sm)
{
    sm.runtime.progress.print();
}

template <typename T>
void progress_clear_line(typename T::StateMachine &sm)
{
    aviize::erase_line(sm.runtime.progress.text);
    sm.runtime.progress.print();
}

template <typename T>
void print_stats(typename T::StateMachine &sm)
{
    aviize::print(sm.runtime.stats());
}

template <typename T>
void create_stats_file(typename T::StateMachine &sm)
{
    if (sm.settings.stats_file_name.empty()) {
        return;
    }

    std::ofstream stats_f(sm.settings.stats_file_name);
    stats_f << sm.runtimeget_stats();
}

template <typename T>
void init_rec_periods(typename T::StateMachine &sm)
{
    // record state queue holds the
    // numbers of all states that must produce a record
    // - always starts with 1
    // - always end n_cycles
    assert(sm.runtime.rec_cycles_queue.empty());
    for (size_t i = 1; i <= sm.settings.n_cycles; i++) {
        if (sm.settings.rec_schedule.is_time(i)) {
            sm.runtime.rec_cycles_queue.push(i);
        }
    }
}

template <typename T>
void decide_rec(typename T::StateMachine &sm)
{
    if (!sm.runtime.rec_cycles_queue.empty() &&
        sm.runtime.cycle_i >= sm.runtime.rec_cycles_queue.front()) {
        sm.runtime.rec_cycles_queue.pop();
        sm.runtime.do_rec = true;
        return;
    }

    sm.runtime.do_rec = false;
}

template <typename T>
void decide_smfg_done(typename T::StateMachine &sm)
{
    assert(!sm.runtime.smfg_done);
    if (sm.runtime.cycle_i == sm.settings.n_cycles) {
        sm.runtime.smfg_done = true;
    }
}

template <typename T>
void cycle_i_inc(typename T::StateMachine &sm)
{
    sm.cycle_i++;
}

}  // namespace opiitis
