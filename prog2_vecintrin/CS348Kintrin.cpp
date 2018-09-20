#include "CS348Kintrin.h"
#include "logger.h"

//******************
//* Implementation *
//******************

__cs348k_mask _cs348k_init_ones(int first) {
  __cs348k_mask mask;
  for (int i=0; i<VECTOR_WIDTH; i++) {
    mask.value[i] = (i<first) ? true : false;
  }
  return mask;
}

__cs348k_mask _cs348k_mask_not(__cs348k_mask &maska) {
  __cs348k_mask resultMask;
  for (int i=0; i<VECTOR_WIDTH; i++) {
    resultMask.value[i] = !maska.value[i];
  }
  CS348KLogger.addLog("masknot", _cs348k_init_ones(), VECTOR_WIDTH);
  return resultMask;
}

__cs348k_mask _cs348k_mask_or(__cs348k_mask &maska, __cs348k_mask &maskb) {
  __cs348k_mask resultMask;
  for (int i=0; i<VECTOR_WIDTH; i++) {
    resultMask.value[i] = maska.value[i] | maskb.value[i];
  }
  CS348KLogger.addLog("maskor", _cs348k_init_ones(), VECTOR_WIDTH);
  return resultMask;
}

__cs348k_mask _cs348k_mask_and(__cs348k_mask &maska, __cs348k_mask &maskb) {
  __cs348k_mask resultMask;
  for (int i=0; i<VECTOR_WIDTH; i++) {
    resultMask.value[i] = maska.value[i] && maskb.value[i];
  }
  CS348KLogger.addLog("maskand", _cs348k_init_ones(), VECTOR_WIDTH);
  return resultMask;
}

int _cs348k_cntbits(__cs348k_mask &maska) {
  int count = 0;
  for (int i=0; i<VECTOR_WIDTH; i++) {
    if (maska.value[i]) count++;
  }
  CS348KLogger.addLog("cntbits", _cs348k_init_ones(), VECTOR_WIDTH);
  return count;
}

template <typename T>
void _cs348k_vset(__cs348k_vec<T> &vecResult, T value, __cs348k_mask &mask) {
  for (int i=0; i<VECTOR_WIDTH; i++) {
    vecResult.value[i] = mask.value[i] ? value : vecResult.value[i];
  }
  CS348KLogger.addLog("vset", mask, VECTOR_WIDTH);
}

template void _cs348k_vset<float>(__cs348k_vec_float &vecResult, float value, __cs348k_mask &mask);
template void _cs348k_vset<int>(__cs348k_vec_int &vecResult, int value, __cs348k_mask &mask);

void _cs348k_vset_float(__cs348k_vec_float &vecResult, float value, __cs348k_mask &mask) { _cs348k_vset<float>(vecResult, value, mask); }
void _cs348k_vset_int(__cs348k_vec_int &vecResult, int value, __cs348k_mask &mask) { _cs348k_vset<int>(vecResult, value, mask); }

__cs348k_vec_float _cs348k_vset_float(float value) {
  __cs348k_vec_float vecResult;
  __cs348k_mask mask = _cs348k_init_ones();
  _cs348k_vset_float(vecResult, value, mask);
  return vecResult;
}
__cs348k_vec_int _cs348k_vset_int(int value) {
  __cs348k_vec_int vecResult;
  __cs348k_mask mask = _cs348k_init_ones();
  _cs348k_vset_int(vecResult, value, mask);
  return vecResult;
}

template <typename T>
void _cs348k_vmove(__cs348k_vec<T> &dest, __cs348k_vec<T> &src, __cs348k_mask &mask) {
    for (int i = 0; i < VECTOR_WIDTH; i++) {
        dest.value[i] = mask.value[i] ? src.value[i] : dest.value[i];
    }
    CS348KLogger.addLog("vmove", mask, VECTOR_WIDTH);
}

template void _cs348k_vmove<float>(__cs348k_vec_float &dest, __cs348k_vec_float &src, __cs348k_mask &mask);
template void _cs348k_vmove<int>(__cs348k_vec_int &dest, __cs348k_vec_int &src, __cs348k_mask &mask);

void _cs348k_vmove_float(__cs348k_vec_float &dest, __cs348k_vec_float &src, __cs348k_mask &mask) { _cs348k_vmove<float>(dest, src, mask); }
void _cs348k_vmove_int(__cs348k_vec_int &dest, __cs348k_vec_int &src, __cs348k_mask &mask) { _cs348k_vmove<int>(dest, src, mask); }

template <typename T>
void _cs348k_vload(__cs348k_vec<T> &dest, T* src, __cs348k_mask &mask) {
  for (int i=0; i<VECTOR_WIDTH; i++) {
    dest.value[i] = mask.value[i] ? src[i] : dest.value[i];
  }
  CS348KLogger.addLog("vload", mask, VECTOR_WIDTH);
}

template void _cs348k_vload<float>(__cs348k_vec_float &dest, float* src, __cs348k_mask &mask);
template void _cs348k_vload<int>(__cs348k_vec_int &dest, int* src, __cs348k_mask &mask);

void _cs348k_vload_float(__cs348k_vec_float &dest, float* src, __cs348k_mask &mask) { _cs348k_vload<float>(dest, src, mask); }
void _cs348k_vload_int(__cs348k_vec_int &dest, int* src, __cs348k_mask &mask) { _cs348k_vload<int>(dest, src, mask); }

template <typename T>
void _cs348k_vstore(T* dest, __cs348k_vec<T> &src, __cs348k_mask &mask) {
  for (int i=0; i<VECTOR_WIDTH; i++) {
    dest[i] = mask.value[i] ? src.value[i] : dest[i];
  }
  CS348KLogger.addLog("vstore", mask, VECTOR_WIDTH);
}

template void _cs348k_vstore<float>(float* dest, __cs348k_vec_float &src, __cs348k_mask &mask);
template void _cs348k_vstore<int>(int* dest, __cs348k_vec_int &src, __cs348k_mask &mask);

void _cs348k_vstore_float(float* dest, __cs348k_vec_float &src, __cs348k_mask &mask) { _cs348k_vstore<float>(dest, src, mask); }
void _cs348k_vstore_int(int* dest, __cs348k_vec_int &src, __cs348k_mask &mask) { _cs348k_vstore<int>(dest, src, mask); }

template <typename T>
void _cs348k_vadd(__cs348k_vec<T> &vecResult, __cs348k_vec<T> &veca, __cs348k_vec<T> &vecb, __cs348k_mask &mask) {
  for (int i=0; i<VECTOR_WIDTH; i++) {
    vecResult.value[i] = mask.value[i] ? (veca.value[i] + vecb.value[i]) : vecResult.value[i];
  }
  CS348KLogger.addLog("vadd", mask, VECTOR_WIDTH);
}

template void _cs348k_vadd<float>(__cs348k_vec_float &vecResult, __cs348k_vec_float &veca, __cs348k_vec_float &vecb, __cs348k_mask &mask);
template void _cs348k_vadd<int>(__cs348k_vec_int &vecResult, __cs348k_vec_int &veca, __cs348k_vec_int &vecb, __cs348k_mask &mask);

void _cs348k_vadd_float(__cs348k_vec_float &vecResult, __cs348k_vec_float &veca, __cs348k_vec_float &vecb, __cs348k_mask &mask) { _cs348k_vadd<float>(vecResult, veca, vecb, mask); }
void _cs348k_vadd_int(__cs348k_vec_int &vecResult, __cs348k_vec_int &veca, __cs348k_vec_int &vecb, __cs348k_mask &mask) { _cs348k_vadd<int>(vecResult, veca, vecb, mask); }

template <typename T>
void _cs348k_vsub(__cs348k_vec<T> &vecResult, __cs348k_vec<T> &veca, __cs348k_vec<T> &vecb, __cs348k_mask &mask) {
  for (int i=0; i<VECTOR_WIDTH; i++) {
    vecResult.value[i] = mask.value[i] ? (veca.value[i] - vecb.value[i]) : vecResult.value[i];
  }
  CS348KLogger.addLog("vsub", mask, VECTOR_WIDTH);
}

template void _cs348k_vsub<float>(__cs348k_vec_float &vecResult, __cs348k_vec_float &veca, __cs348k_vec_float &vecb, __cs348k_mask &mask);
template void _cs348k_vsub<int>(__cs348k_vec_int &vecResult, __cs348k_vec_int &veca, __cs348k_vec_int &vecb, __cs348k_mask &mask);

void _cs348k_vsub_float(__cs348k_vec_float &vecResult, __cs348k_vec_float &veca, __cs348k_vec_float &vecb, __cs348k_mask &mask) { _cs348k_vsub<float>(vecResult, veca, vecb, mask); }
void _cs348k_vsub_int(__cs348k_vec_int &vecResult, __cs348k_vec_int &veca, __cs348k_vec_int &vecb, __cs348k_mask &mask) { _cs348k_vsub<int>(vecResult, veca, vecb, mask); }

template <typename T>
void _cs348k_vmult(__cs348k_vec<T> &vecResult, __cs348k_vec<T> &veca, __cs348k_vec<T> &vecb, __cs348k_mask &mask) {
  for (int i=0; i<VECTOR_WIDTH; i++) {
    vecResult.value[i] = mask.value[i] ? (veca.value[i] * vecb.value[i]) : vecResult.value[i];
  }
  CS348KLogger.addLog("vmult", mask, VECTOR_WIDTH);
}

template void _cs348k_vmult<float>(__cs348k_vec_float &vecResult, __cs348k_vec_float &veca, __cs348k_vec_float &vecb, __cs348k_mask &mask);
template void _cs348k_vmult<int>(__cs348k_vec_int &vecResult, __cs348k_vec_int &veca, __cs348k_vec_int &vecb, __cs348k_mask &mask);

void _cs348k_vmult_float(__cs348k_vec_float &vecResult, __cs348k_vec_float &veca, __cs348k_vec_float &vecb, __cs348k_mask &mask) { _cs348k_vmult<float>(vecResult, veca, vecb, mask); }
void _cs348k_vmult_int(__cs348k_vec_int &vecResult, __cs348k_vec_int &veca, __cs348k_vec_int &vecb, __cs348k_mask &mask) { _cs348k_vmult<int>(vecResult, veca, vecb, mask); }

template <typename T>
void _cs348k_vdiv(__cs348k_vec<T> &vecResult, __cs348k_vec<T> &veca, __cs348k_vec<T> &vecb, __cs348k_mask &mask) {
  for (int i=0; i<VECTOR_WIDTH; i++) {
    vecResult.value[i] = mask.value[i] ? (veca.value[i] / vecb.value[i]) : vecResult.value[i];
  }
  CS348KLogger.addLog("vdiv", mask, VECTOR_WIDTH);
}

template void _cs348k_vdiv<float>(__cs348k_vec_float &vecResult, __cs348k_vec_float &veca, __cs348k_vec_float &vecb, __cs348k_mask &mask);
template void _cs348k_vdiv<int>(__cs348k_vec_int &vecResult, __cs348k_vec_int &veca, __cs348k_vec_int &vecb, __cs348k_mask &mask);

void _cs348k_vdiv_float(__cs348k_vec_float &vecResult, __cs348k_vec_float &veca, __cs348k_vec_float &vecb, __cs348k_mask &mask) { _cs348k_vdiv<float>(vecResult, veca, vecb, mask); }
void _cs348k_vdiv_int(__cs348k_vec_int &vecResult, __cs348k_vec_int &veca, __cs348k_vec_int &vecb, __cs348k_mask &mask) { _cs348k_vdiv<int>(vecResult, veca, vecb, mask); }

template <typename T>
void _cs348k_vabs(__cs348k_vec<T> &vecResult, __cs348k_vec<T> &veca, __cs348k_mask &mask) {
  for (int i=0; i<VECTOR_WIDTH; i++) {
    vecResult.value[i] = mask.value[i] ? (abs(veca.value[i])) : vecResult.value[i];
  }
  CS348KLogger.addLog("vabs", mask, VECTOR_WIDTH);
}

template void _cs348k_vabs<float>(__cs348k_vec_float &vecResult, __cs348k_vec_float &veca, __cs348k_mask &mask);
template void _cs348k_vabs<int>(__cs348k_vec_int &vecResult, __cs348k_vec_int &veca, __cs348k_mask &mask);

void _cs348k_vabs_float(__cs348k_vec_float &vecResult, __cs348k_vec_float &veca, __cs348k_mask &mask) { _cs348k_vabs<float>(vecResult, veca, mask); }
void _cs348k_vabs_int(__cs348k_vec_int &vecResult, __cs348k_vec_int &veca, __cs348k_mask &mask) { _cs348k_vabs<int>(vecResult, veca, mask); }

template <typename T>
void _cs348k_vgt(__cs348k_mask &maskResult, __cs348k_vec<T> &veca, __cs348k_vec<T> &vecb, __cs348k_mask &mask) {
  for (int i=0; i<VECTOR_WIDTH; i++) {
    maskResult.value[i] = mask.value[i] ? (veca.value[i] > vecb.value[i]) : maskResult.value[i];
  }
  CS348KLogger.addLog("vgt", mask, VECTOR_WIDTH);
}

template void _cs348k_vgt<float>(__cs348k_mask &maskResult, __cs348k_vec_float &veca, __cs348k_vec_float &vecb, __cs348k_mask &mask);
template void _cs348k_vgt<int>(__cs348k_mask &maskResult, __cs348k_vec_int &veca, __cs348k_vec_int &vecb, __cs348k_mask &mask);

void _cs348k_vgt_float(__cs348k_mask &maskResult, __cs348k_vec_float &veca, __cs348k_vec_float &vecb, __cs348k_mask &mask) { _cs348k_vgt<float>(maskResult, veca, vecb, mask); }
void _cs348k_vgt_int(__cs348k_mask &maskResult, __cs348k_vec_int &veca, __cs348k_vec_int &vecb, __cs348k_mask &mask) { _cs348k_vgt<int>(maskResult, veca, vecb, mask); }

template <typename T>
void _cs348k_vlt(__cs348k_mask &maskResult, __cs348k_vec<T> &veca, __cs348k_vec<T> &vecb, __cs348k_mask &mask) {
  for (int i=0; i<VECTOR_WIDTH; i++) {
    maskResult.value[i] = mask.value[i] ? (veca.value[i] < vecb.value[i]) : maskResult.value[i];
  }
  CS348KLogger.addLog("vlt", mask, VECTOR_WIDTH);
}

template void _cs348k_vlt<float>(__cs348k_mask &maskResult, __cs348k_vec_float &veca, __cs348k_vec_float &vecb, __cs348k_mask &mask);
template void _cs348k_vlt<int>(__cs348k_mask &maskResult, __cs348k_vec_int &veca, __cs348k_vec_int &vecb, __cs348k_mask &mask);

void _cs348k_vlt_float(__cs348k_mask &maskResult, __cs348k_vec_float &veca, __cs348k_vec_float &vecb, __cs348k_mask &mask) { _cs348k_vlt<float>(maskResult, veca, vecb, mask); }
void _cs348k_vlt_int(__cs348k_mask &maskResult, __cs348k_vec_int &veca, __cs348k_vec_int &vecb, __cs348k_mask &mask) { _cs348k_vlt<int>(maskResult, veca, vecb, mask); }

template <typename T>
void _cs348k_veq(__cs348k_mask &maskResult, __cs348k_vec<T> &veca, __cs348k_vec<T> &vecb, __cs348k_mask &mask) {
  for (int i=0; i<VECTOR_WIDTH; i++) {
    maskResult.value[i] = mask.value[i] ? (veca.value[i] == vecb.value[i]) : maskResult.value[i];
  }
  CS348KLogger.addLog("veq", mask, VECTOR_WIDTH);
}

template void _cs348k_veq<float>(__cs348k_mask &maskResult, __cs348k_vec_float &veca, __cs348k_vec_float &vecb, __cs348k_mask &mask);
template void _cs348k_veq<int>(__cs348k_mask &maskResult, __cs348k_vec_int &veca, __cs348k_vec_int &vecb, __cs348k_mask &mask);

void _cs348k_veq_float(__cs348k_mask &maskResult, __cs348k_vec_float &veca, __cs348k_vec_float &vecb, __cs348k_mask &mask) { _cs348k_veq<float>(maskResult, veca, vecb, mask); }
void _cs348k_veq_int(__cs348k_mask &maskResult, __cs348k_vec_int &veca, __cs348k_vec_int &vecb, __cs348k_mask &mask) { _cs348k_veq<int>(maskResult, veca, vecb, mask); }

template <typename T>
void _cs348k_hadd(__cs348k_vec<T> &vecResult, __cs348k_vec<T> &vec) {
  for (int i=0; i<VECTOR_WIDTH/2; i++) {
    T result = vec.value[2*i] + vec.value[2*i+1];
    vecResult.value[2 * i] = result;
    vecResult.value[2 * i + 1] = result;
  }
}

template void _cs348k_hadd<float>(__cs348k_vec_float &vecResult, __cs348k_vec_float &vec);

void _cs348k_hadd_float(__cs348k_vec_float &vecResult, __cs348k_vec_float &vec) { _cs348k_hadd<float>(vecResult, vec); }

template <typename T>
void _cs348k_interleave(__cs348k_vec<T> &vecResult, __cs348k_vec<T> &vec) {
  for (int i=0; i<VECTOR_WIDTH; i++) {
    int index = i < VECTOR_WIDTH/2 ? (2 * i) : (2 * (i - VECTOR_WIDTH/2) + 1);
    vecResult.value[i] = vec.value[index];
  }
}

template void _cs348k_interleave<float>(__cs348k_vec_float &vecResult, __cs348k_vec_float &vec);

void _cs348k_interleave_float(__cs348k_vec_float &vecResult, __cs348k_vec_float &vec) { _cs348k_interleave<float>(vecResult, vec); }

void addUserLog(const char * logStr) {
  CS348KLogger.addLog(logStr, _cs348k_init_ones(), 0);
}

