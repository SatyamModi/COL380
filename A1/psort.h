#pragma once

#include <stdint.h>

void SequentialSort(uint32_t *data, uint32_t n);
void ParallelSort(uint32_t *data, uint32_t n, int p);