#ifdef __cplusplus
extern "C" {
#endif

/*
 * Copyright (C) 2013 Mark Hills <mark@xwax.org>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * version 2, as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License version 2 for more details.
 *
 * You should have received a copy of the GNU General Public License
 * version 2 along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 *
 */
// ref :http://www.pogo.org.uk/~mark/bpm-tools/
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

//ref: https://github.com/mackron/dr_libs/blob/master/dr_wav.h
#define DR_WAV_IMPLEMENTATION

#include "dr_wav.h"

float *wavRead_f32(char *filename, uint32_t *sampleRate, uint64_t *totalSampleCount, unsigned int *channels) {
    float *buffer = drwav_open_and_read_file_f32(filename, channels, sampleRate, totalSampleCount);
    if (buffer == NULL) {
        fprintf(stderr, "read file error.\n");
        exit(1);
    }
    if (*channels == 2) {
        float *bufferSave = buffer;
        for (uint64_t i = 0; i < *totalSampleCount; i += 2) {
            *bufferSave++ = ((buffer[i] + buffer[i + 1]) * 0.5f);
        }
        *totalSampleCount = *totalSampleCount >> 1;
        *channels = 1;
    } else if (*channels != 1) {
        drwav_free(buffer);
        buffer = NULL;
        *sampleRate = 0;
        *totalSampleCount = 0;
    }
    return buffer;
}


#ifndef HAVEFP
#define HAVEFP 1
#define MASK    ((unsigned)(1 << (16 - 1)) + (1 << (16 - 1)) - 1)
#define LOW(x)    ((unsigned)(x) & MASK)
#define HIGH(x)    LOW((x) >> 16)
#define MUL(x, y, z)    { long l = (long)(x) * (long)(y); \
        (z)[0] = LOW(l); (z)[1] = HIGH(l); }
#define CARRY(x, y)    ((long)(x) + (long)(y) > MASK)
#define ADDEQU(x, y, z)    ((z) = CARRY(x, (y)), (x) = LOW((x) + (y)))
static unsigned x[3] = {0x330E, 0xABCD, 0x1234}, a[3] = {0xE66D, 0xDEEC, 0x5}, c = 0xB;

static void next() {
    unsigned p[2], q[2], r[2], carry0, carry1;
    MUL(a[0], x[0], p);
    ADDEQU(p[0], c, carry0);
    ADDEQU(p[1], carry0, carry1);
    MUL(a[0], x[1], q);
    ADDEQU(p[1], q[0], carry0);
    MUL(a[1], x[0], r);
    x[2] = LOW(carry0 + carry1 + CARRY(p[1], r[0]) + q[1] + r[1] +
               a[0] * x[2] + a[1] * x[1] + a[2] * x[0]);
    x[1] = LOW(p[1] + r[0]);
    x[0] = LOW(p[0]);
}

double
drand48() {
    static double two16m = 1.0 / (1L << 16);
    next();
    return (two16m * (two16m * (two16m * x[0] + x[1]) + x[2]));
}

#endif


/*
 * Sample from the metered energy
 *
 * No need to interpolate and it makes a tiny amount of difference; we
 * take a random sample of samples, any errors are averaged out.
 */

static double sample(float nrg[], size_t len, double offset) {
    double n;
    size_t i;

    n = floor(offset);
    i = (size_t) n;

    return (n >= 0.0 && n < (double) len) ? nrg[i] : 0.0;
}

/*
 * Test an autodifference for the given interval
 */

double autodifference(float nrg[], size_t len, double interval) {
    size_t n;
    double mid, v, diff, total;
    static const double beats[] = {-32, -16, -8, -4, -2, -1,
                                   1, 2, 4, 8, 16, 32},
            nobeats[] = {-0.5, -0.25, 0.25, 0.5};

    mid = drand48() * len;
    v = sample(nrg, len, mid);

    diff = 0.0;
    total = 0.0;

    for (n = 0; n < (sizeof(beats) / sizeof(*(beats))); n++) {
        double y, w;

        y = sample(nrg, len, mid + beats[n] * interval);

        w = 1.0 / fabs(beats[n]);
        diff += w * fabs(y - v);
        total += w;
    }

    for (n = 0; n < (sizeof(nobeats) / sizeof(*(nobeats))); n++) {
        double y, w;

        y = sample(nrg, len, mid + nobeats[n] * interval);

        w = fabs(nobeats[n]);
        diff -= w * fabs(y - v);
        total += w;
    }

    return diff / total;
}

/*
 * Beats-per-minute to a sampling interval in energy space
 */

double bpm_to_interval(double bpm, int sampleRate, int interval) {
    double beats_per_second, samples_per_beat;

    beats_per_second = bpm / 60;
    samples_per_beat = sampleRate / beats_per_second;
    return samples_per_beat / interval;
}

/*
 * Sampling interval in enery space to beats-per-minute
 */

double interval_to_bpm(double cur_interval, int interval, int sampleRate) {
    double samples_per_beat, beats_per_second;

    samples_per_beat = cur_interval * interval;
    beats_per_second = (double) sampleRate / samples_per_beat;
    return beats_per_second * 60;
}

/*
 * Scan a range of BPM values for the one with the
 * minimum autodifference
 */
double scan_for_bpm(float nrg[], size_t len,
                    double slowest, double fastest,
                    unsigned int steps,
                    unsigned int samples,
                    int sampleRate, int interval) {
    double step, trough, height;
    unsigned int s;

    slowest = bpm_to_interval(slowest, sampleRate, interval);
    fastest = bpm_to_interval(fastest, sampleRate, interval);
    step = (slowest - fastest) / steps;

    height = INFINITY;
    trough = NAN;

    for (double cur_interval = fastest; cur_interval <= slowest; cur_interval += step) {
        double t = 0.0;
        for (s = 0; s < samples; s++)
            t += autodifference(nrg, len, cur_interval);

        printf("%lf\t%lf\n", interval_to_bpm(cur_interval, interval, sampleRate), t / samples);

        /* Track the lowest value */
        if (t < height) {
            trough = cur_interval;
            height = t;
        }
    }
    return interval_to_bpm(trough, interval, sampleRate);
}

void bpm(const float *data_in, uint32_t sampleRate, uint64_t totalSampleCount, double min, double max) {
    int interval = 128;
    size_t len = totalSampleCount / interval;
    float *nrg = (float *) calloc(len, sizeof(float));
    if (nrg) {

        float v = 0.0;
        int n = 0;
        int index = 0;
        for (int i = 0; i < totalSampleCount; ++i) {
            /* Maintain an energy meter (similar to PPM) */
            float z = data_in[i];
            z = fabsf(z);
            if (z > v) {
                v += (z - v) / 8;
            } else {
                v -= (v - z) / 512;
            }
            /* At regular intervals, sample the energy to give a
             * low-resolution overview of the track */
            n++;
            if (n != interval)
                continue;
            nrg[index++] = v;
            n = 0;
        }
        double bpm = scan_for_bpm(nrg, len, min, max, 1024, 1024, sampleRate, interval);
        printf("bpm: %0.3f \n", bpm);
        free(nrg);
    }
}

int main(int argc, char *argv[]) {
    printf("Audio Processing\n");
    printf("blog:http://cpuimage.cnblogs.com/\n");
    printf("tempo analysis (bpm detection) algorithms in C\n");
    if (argc < 2)
        return -1;

    char *in_file = argv[1];
    uint32_t sampleRate = 0;
    uint64_t totalSampleCount = 0;
    uint32_t channels = 0;
    double min = 84.0, max = 146.0;
    float *data_in = wavRead_f32(in_file, &sampleRate, &totalSampleCount, &channels);
    if (data_in) {
        bpm(data_in, sampleRate, totalSampleCount, min, max);
        free(data_in);
    }
    printf("press any key to exit.\n");
    getchar();
    return 0;
}

#ifdef __cplusplus
}
#endif