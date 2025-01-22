#ifndef FIDITI_UTIL_TIME_HPP
#define FIDITI_UTIL_TIME_HPP

#include <chrono>
#include <ctime>
#include <unordered_map>
#include <string>


namespace fidi
{

class Timings
{
    using RtClock = std::chrono::high_resolution_clock;

    struct Time
    {
        // Milliseconds: realTime.count()
        std::chrono::duration<double, std::milli> realTime;

        // Milliseconds: cpuTime * 1000.0 / CLOCKS_PER_SEC
        std::clock_t cpuTime = 0;
    };

    std::unordered_map<std::string, Time> times;

    std::chrono::time_point<RtClock> rt;
    std::clock_t ct = 0;

 public:
    // Restart the clock (before a sequence of operations to measure)
    void start()
    {
        rt = RtClock::now();
        ct = std::clock();
    }

    // Add elapsed time to a particular counter, and reset the clock
    void add(std::string id)
    {
        auto rt_now = RtClock::now();
        auto ct_now = std::clock();
        times[id].realTime += rt_now - rt;
        times[id].cpuTime += ct_now - ct;

        start();
    }

    void reset()
    {
        times.clear();
        start();
    }

    Time operator[](std::string id) const
    {
        auto it = times.find(id);
        assert(it != times.end());
        return it->second;
    }

    double real_ms(std::string id) const
    {
        return (*this)[id].realTime.count();
    }

    double cpu_ms(std::string id) const
    {
        return (*this)[id].cpuTime * 1000.0 / CLOCKS_PER_SEC;
    }
};


} // namespace fidi


#endif // FIDITI_UTIL_TIME_HPP
