#ifndef SCHISM_TIMING_MACROS_INC
#define SCHISM_TIMING_MACROS_INC

#ifdef INCLUDE_TIMING
#define SCHISM_TIMING_START(started_at) call schism_timing_start_impl(started_at)
#define SCHISM_TIMING_ACCUMULATE(accumulator, started_at) call schism_timing_accumulate_impl(accumulator, started_at)
#else
#define SCHISM_TIMING_START(started_at)
#define SCHISM_TIMING_ACCUMULATE(accumulator, started_at)
#endif

#ifdef TIMER2
#define SCHISM_TIMER2_START(started_at) call schism_timing_start_impl(started_at)
#define SCHISM_TIMER2_CHECKPOINT(label, step, started_at) call schism_timer2_checkpoint_impl(label, step, started_at)
#else
#define SCHISM_TIMER2_START(started_at)
#define SCHISM_TIMER2_CHECKPOINT(label, step, started_at)
#endif

#endif
