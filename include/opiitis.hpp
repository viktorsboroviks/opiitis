#include <boost/json/src.hpp>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include "aviize.hpp"
#include "garaaza.hpp"
#include "iestaade.hpp"

namespace opiitis {

// TODO: review
// forward declarations
// clang-format off
struct Settings;
template <typename T> struct StateParams;
template <typename T> struct State;
template <typename T> struct StateMachine;
// clang-format on

// TODO: review
// all classes packaged as a single template parameter
struct Context {
    using Settings    = Settings;
    using StateParams = StateParams<Context>;
    using State       = State<Context>;
};

// TODO: move somewhere else
// tools
bool last_line_empty(const std::string &text)
{
    return (text.empty() || text.back() == '\n' || text.back() == '\r');
}

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
    size_t n_states = 1000;

    // pools
    size_t superpool_size  = 100;
    size_t pool_size_min   = 10;
    double pool_inc_coef   = 2.0;
    double pool_kstest_min = 0.01;

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
        n_states                = iestaade::from_json<size_t>     (config_path, key_prefix + "n_states");

        superpool_size          = iestaade::from_json<size_t>     (config_path, key_prefix + "superpool_size");
        pool_size_min           = iestaade::from_json<size_t>     (config_path, key_prefix + "pool_size_min");
        pool_inc_coef           = iestaade::from_json<double>     (config_path, key_prefix + "pool_inc_coef");
        pool_kstest_min         = iestaade::from_json<double>     (config_path, key_prefix + "pool_kstest_min");

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
    // run status
    bool do_rec    = false;
    size_t cycle_i = 1;
    size_t smfg_i  = 0;
    bool smfg_done = false;

    // TODO: do we need it for the base class?
    // state
    typename T::State state;

    // performance
    std::chrono::time_point<std::chrono::steady_clock> start_time;
    double cycle_time_us = 0;

    // history
    std::queue<size_t> rec_states_queue;

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
};

// TODO: review
template <typename T = Context>
struct StateParams {
    virtual void randomize() {}

    virtual void change() {}
};

// TODO: review
template <typename T = Context>
struct State {
    State()
    {
        // TODO: consider removal
    }

    explicit State(const typename T::StateParams &params)
    {
        (void)params;
        // TODO: consider removal
    }

    virtual double get_energy()
    {
        return -1;
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
        std::string name;
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

// TODO: review all below
template <typename T>
void init_log(StateMachine<T> &sm)
{
    assert(!sm.log_f.is_open());
    if (sm.settings.log_file_name.empty()) {
        return;
    }

    sm.log_f.open(sm.settings.log_file_name);
    sm.log_f << "run_i,";
    sm.log_f << "t,";
    sm.log_f << "e";
    sm.log_f << std::endl;
}

template <typename T>
void update_log(StateMachine<T> &sm)
{
    // write the log to file
    assert(!sm.settings.log_file_name.empty());
    if (!sm.log_f.is_open()) {
        return;
    }

    sm.log_f << sm.run_i;
    sm.log_f << "," << sm.state.get_energy();
    sm.log_f << std::endl;
}

template <typename T>
void progress_init_init_loop(StateMachine<T> &sm)
{
    sm.progress.n_min         = 1;
    sm.progress.n_max         = sm.settings.init_t_max_attempts;
    sm.progress.update_period = sm.settings.progress_update_period;
}

template <typename T>
void progress_init_run_loop(StateMachine<T> &sm)
{
    sm.progress.n_min         = 1;
    sm.progress.n_max         = sm.settings.n_states;
    sm.progress.update_period = sm.settings.progress_update_period;
}

template <typename T>
void progress_text_reset(StateMachine<T> &sm)
{
    sm.progress.reset();
}

template <typename T>
void progress_text_add_nl(StateMachine<T> &sm)
{
    sm.progress.text += "\n";
}

template <typename T>
void progress_text_add_stats(StateMachine<T> &sm)
{
    const double run_s = 1 / sm.cycle_time_us * 1000000;

    std::stringstream ss;
    ss << " n/s " << run_s;

    sm.progress.text += std::string(ss.str());
}

template <typename T>
void progress_text_add_total(StateMachine<T> &sm)
{
    if (!last_line_empty(sm.progress.text)) {
        sm.progress.text += " ";
    }
    sm.progress.text += sm.progress.str_total(sm.run_i);
}

template <typename T>
void progress_text_add_pct(StateMachine<T> &sm)
{
    if (!last_line_empty(sm.progress.text)) {
        sm.progress.text += " ";
    }
    sm.progress.text += sm.progress.str_pct(sm.run_i);
}

template <typename T>
void progress_text_add_eta(StateMachine<T> &sm)
{
    if (!last_line_empty(sm.progress.text)) {
        sm.progress.text += " ";
    }
    sm.progress.text += sm.progress.str_eta(sm.run_i);
}

template <typename T>
void progress_text_add_energy(StateMachine<T> &sm)
{
    if (!last_line_empty(sm.progress.text)) {
        sm.progress.text += " ";
    }
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3) << sm.state.get_energy();
    sm.progress.text += "e " + oss.str();
}

template <typename T>
void progress_print(StateMachine<T> &sm)
{
    sm.progress.print();
}

template <typename T>
void progress_clear_line(StateMachine<T> &sm)
{
    aviize::erase_line(sm.progress.text);
    sm.progress.print();
}

template <typename T>
void print_stats(StateMachine<T> &sm)
{
    std::stringstream ss{};
    ss << sm.get_stats();
    aviize::print(ss);
}

template <typename T>
void create_stats_file(StateMachine<T> &sm)
{
    if (sm.settings.stats_file_name.empty()) {
        return;
    }

    std::ofstream f(sm.settings.stats_file_name);
    f << sm.get_stats();
}

template <typename T>
void init_rec_periods(StateMachine<T> &sm)
{
    // record state queue holds the
    // numbers of all states that must produce a record
    // - always starts with 1
    // - always end n_states
    assert(sm.rec_states_queue.empty());
    for (size_t i = 1; i <= sm.settings.n_states; i++) {
        if (sm.settings.rec_schedule.is_time(i)) {
            sm.rec_states_queue.push(i);
        }
    }
}

template <typename T>
void decide_rec(StateMachine<T> &sm)
{
    if (!sm.rec_states_queue.empty() &&
        sm.run_i >= sm.rec_states_queue.front()) {
        sm.rec_states_queue.pop();
        sm.do_rec = true;
        return;
    }

    sm.do_rec = false;
}

template <typename T>
void decide_run_done(StateMachine<T> &sm)
{
    assert(!sm.run_done);
    if (sm.run_i == sm.settings.n_states) {
        sm.run_done = true;
    }
}

template <typename T>
void run_i_inc(StateMachine<T> &sm)
{
    sm.run_i++;
}

}  // namespace opiitis
