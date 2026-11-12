#ifndef GPUTIMER_H
#define GPUTIMER_H

// GPU time per render pass, from GL timestamp queries (benchmarks only: --bench prints it).
// GpuTimer_Mark(name) stamps the point where the named pass ENDS; the time charged to it is from
// the previous mark in the same frame. Results are read back a few frames later, so they never
// stall the pipeline. Off (no GL calls) until GpuTimer_Enable(true).

void GpuTimer_Enable(bool on);
void GpuTimer_BeginFrame();            // before the first draw of a frame
void GpuTimer_Mark(const char* name);  // a string literal: the pass that just finished
void GpuTimer_Print();                 // "bench: gpu ..." averages per pass, then the total
// CPU time on the main thread between marks: GpuTimer_Mark also marks the CPU side under its name.
void GpuTimer_CpuFrameStart();          // top of the frame loop
void GpuTimer_CpuMark(const char* name);
void GpuTimer_PrintCpu();

#endif
