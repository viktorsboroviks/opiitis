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
template <typename T> struct Engine;
// clang-format on

// all classes packaged as a single template parameter
struct Context {
    using Settings = Settings;
    using Runtime  = Runtime<Context>;
    using Engine   = Engine<Context>;
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
    // settings
    const typename T::Settings settings;

    // status
    size_t engine_state_i  = 0;
    bool engine_state_done = false;
    size_t cycle_i         = 0;

    // performance
    std::chrono::time_point<std::chrono::steady_clock> start_time;
    double cycle_time_us = 0;

    // recording
    std::queue<size_t> rec_cycles_queue;
    bool do_rec = false;

    // data display and logging
    aviize::Progress progress;
    std::ofstream log_f;

    Runtime(const std::string &config_path, const std::string &key_prefix) :
        settings(config_path, key_prefix)
    {
    }

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
struct Engine {
    typename T::Runtime runtime;

    typedef std::function<void(typename T::Runtime &)> engine_fn_t;
    struct EngineState {
        std::string name;  // cppcheck-suppress unusedStructMember
        bool repeat_until_done = false;
        std::vector<engine_fn_t> functions;
    };

    // operate via add_state() and add_state_cyclic()
    std::vector<EngineState> engine_states;

    explicit Engine(const std::string &config_path,
                    const std::string &key_prefix) :
        runtime(config_path, key_prefix)
    {
    }

    void add_state(const std::string &name,
                   const std::vector<engine_fn_t> &functions)
    {
        engine_states.push_back(EngineState{name, false, functions});
    }

    void add_state_cyclic(const std::string &name,
                          const std::vector<engine_fn_t> &functions)
    {
        engine_states.push_back(EngineState{name, true, functions});
    }

    void run()
    {
        runtime.start_time = std::chrono::steady_clock::now();

        assert(!engine_states.empty());

        // go over all engine states
        // - each es can be run multiple times until engine_state_done
        // - track engine_state_i
        // - track cycle_time_us to know how long it took to run one cycle
        //   of the engine state
        for (runtime.engine_state_i = 0;
             runtime.engine_state_i < engine_states.size();
             runtime.engine_state_i++) {
            const auto &engine_state = engine_states[runtime.engine_state_i];
            auto cycle_begin_time    = std::chrono::steady_clock::now();

            runtime.engine_state_done = false;
            do {
                // go over all functions in the group
                // - abort on engine_state_done
                for (const engine_fn_t &f : engine_state.functions) {
                    if (runtime.engine_state_done) {
                        break;
                    }
                    f(this->runtime);
                }

                // update the cycle time
                const auto cycle_end_time = std::chrono::steady_clock::now();
                runtime.cycle_time_us =
                        std::chrono::duration_cast<std::chrono::microseconds>(
                                cycle_end_time - cycle_begin_time)
                                .count();
                cycle_begin_time = cycle_end_time;

                // if repeat_until_done is not set, abort
                if (!engine_state.repeat_until_done) {
                    break;
                }
            } while (!runtime.engine_state_done);
        }
    }
};

template <typename T>
void init_log(typename T::Runtime &r)
{
    assert(!r.log_f.is_open());
    if (r.settings.log_file_name.empty()) {
        return;
    }

    r.log_f.open(r.settings.log_file_name);
    r.log_f << r.log_header();
    r.log_f << std::endl;
}

template <typename T>
void update_log(typename T::Runtime &r)
{
    assert(!r.settings.log_file_name.empty());
    if (!r.log_f.is_open()) {
        return;
    }

    r.log_f << r.log_line();
    r.log_f << std::endl;
}

template <typename T>
void progress_init(typename T::Runtime &r)
{
    r.progress.n_min         = 1;
    r.progress.n_max         = r.settings.n_cycles;
    r.progress.update_period = r.settings.progress_update_period;
}

template <typename T>
void progress_text_reset(typename T::Runtime &r)
{
    r.progress.reset();
}

template <typename T>
void progress_text_add_nl(typename T::Runtime &r)
{
    r.progress.text += "\n";
}

template <typename T>
void progress_text_add_speed(typename T::Runtime &r)
{
    const double cycle_s = 1 / r.cycle_time_us * 1000000;

    std::stringstream ss;
    ss << " n/s " << cycle_s;

    r.progress.text += std::string(ss.str());
}

template <typename T>
void progress_text_add_total(typename T::Runtime &r)
{
    if (!r.progress.last_line_empty()) {
        r.progress.text += " ";
    }
    r.progress.text += r.progress.str_total(r.cycle_i);
}

template <typename T>
void progress_text_add_pct(typename T::Runtime &r)
{
    if (!r.progress.last_line_empty()) {
        r.progress.text += " ";
    }
    r.progress.text += r.progress.str_pct(r.cycle_i);
}

template <typename T>
void progress_text_add_eta(typename T::Runtime &r)
{
    if (!r.progress.last_line_empty()) {
        r.progress.text += " ";
    }
    r.progress.text += r.progress.str_eta(r.cycle_i);
}

template <typename T>
void progress_print(typename T::Runtime &r)
{
    r.progress.print();
}

template <typename T>
void progress_clear_line(typename T::Runtime &r)
{
    aviize::erase_line(r.progress.text);
    r.progress.print();
}

template <typename T>
void print_stats(typename T::Runtime &r)
{
    aviize::print(r.stats());
}

template <typename T>
void create_stats_file(typename T::Runtime &r)
{
    if (r.settings.stats_file_name.empty()) {
        return;
    }

    std::ofstream stats_f(r.settings.stats_file_name);
    stats_f << r.stats();
}

template <typename T>
void init_rec_periods(typename T::Runtime &r)
{
    // record state queue holds the
    // numbers of all states that must produce a record
    // - always starts with 1
    // - always end n_cycles
    assert(r.rec_cycles_queue.empty());
    for (size_t i = 1; i <= r.settings.n_cycles; i++) {
        if (r.settings.rec_schedule.is_time(i)) {
            r.rec_cycles_queue.push(i);
        }
    }
}

template <typename T>
void decide_rec(typename T::Runtime &r)
{
    if (!r.rec_cycles_queue.empty() &&
        r.cycle_i >= r.rec_cycles_queue.front()) {
        r.rec_cycles_queue.pop();
        r.do_rec = true;
        return;
    }

    r.do_rec = false;
}

template <typename T>
void decide_done(typename T::Runtime &r)
{
    assert(!r.engine_state_done);
    if (r.cycle_i == r.settings.n_cycles) {
        r.engine_state_done = true;
    }
}

template <typename T>
void cycle_i_inc(typename T::Runtime &r)
{
    r.cycle_i++;
}

}  // namespace opiitis
