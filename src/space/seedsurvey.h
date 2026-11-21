#ifndef SEEDSURVEY_H
#define SEEDSURVEY_H

// --seed-survey COUNT FIRSTSEED OUTDIR [THREADS] [SUBSET] (#532): compares the New World previews of many world
// seeds (land masks, height correlation, with and without longitude shifts / flips). No window. Returns the exit code.
int SeedSurvey_Run(int argc, char** argv);

#endif
