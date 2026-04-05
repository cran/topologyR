// ======================================================================
// bitset_ops.h — shared multi-word bitset infrastructure for topologyR
//
// Provides templatized (compile-time W) and runtime-W bitset operations,
// flat hash sets, flat vectors, and decode helpers used across the
// topology engine and Alexandrov engine.
//
// All functions are inline or template to be safe in a shared header.
// ======================================================================

#ifndef TOPOLOGYR_BITSET_OPS_H_
#define TOPOLOGYR_BITSET_OPS_H_

#include <Rcpp.h>
#include <vector>
#include <cstdint>
#include <cstring>


// ======================================================================
// Section 1: Templatized bitset operations
//
// When W is a compile-time constant (especially W=1), the compiler
// eliminates all loops and branch overhead entirely.
// ======================================================================

template<int W>
struct BitOps {
  static inline void bor(uint64_t* dst, const uint64_t* a,
                         const uint64_t* b) {
    for (int i = 0; i < W; i++) dst[i] = a[i] | b[i];
  }
  static inline void band(uint64_t* dst, const uint64_t* a,
                          const uint64_t* b) {
    for (int i = 0; i < W; i++) dst[i] = a[i] & b[i];
  }
  static inline void bnot(uint64_t* dst, const uint64_t* a,
                          const uint64_t* mask) {
    for (int i = 0; i < W; i++) dst[i] = ~a[i] & mask[i];
  }
  static inline bool equal(const uint64_t* a, const uint64_t* b) {
    for (int i = 0; i < W; i++) if (a[i] != b[i]) return false;
    return true;
  }
  static inline bool is_zero(const uint64_t* a) {
    for (int i = 0; i < W; i++) if (a[i] != 0) return false;
    return true;
  }
  static inline bool is_subset(const uint64_t* a, const uint64_t* b) {
    for (int i = 0; i < W; i++) if ((a[i] & b[i]) != a[i]) return false;
    return true;
  }
  static inline size_t hash(const uint64_t* a) {
    size_t h = 0;
    for (int i = 0; i < W; i++) {
      uint64_t x = a[i];
      x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
      x ^= x >> 27; x *= 0x94d049bb133111ebULL;
      x ^= x >> 31;
      h ^= static_cast<size_t>(x) + 0x9e3779b9 + (h << 6) + (h >> 2);
    }
    return h;
  }
  static inline void copy(uint64_t* dst, const uint64_t* src) {
    for (int i = 0; i < W; i++) dst[i] = src[i];
  }
  static inline void zero(uint64_t* dst) {
    for (int i = 0; i < W; i++) dst[i] = 0;
  }
};

inline bool bs_test_bit(const uint64_t* a, int bit) {
  return (a[bit / 64] >> (bit % 64)) & 1ULL;
}

inline void bs_set_bit(uint64_t* a, int bit) {
  a[bit / 64] |= (1ULL << (bit % 64));
}


// ======================================================================
// Section 2: Templatized FlatBitsetVec
// ======================================================================

template<int W>
struct FlatVec {
  std::vector<uint64_t> data;

  int size() const { return static_cast<int>(data.size()) / W; }

  uint64_t* get(int i) { return &data[static_cast<size_t>(i) * W]; }
  const uint64_t* get(int i) const { return &data[static_cast<size_t>(i) * W]; }

  void push(const uint64_t* bs) { data.insert(data.end(), bs, bs + W); }
  void clear() { data.clear(); }
  void reserve(int count) { data.reserve(static_cast<size_t>(count) * W); }
};


// ======================================================================
// Section 3: Templatized FlatHashSet
// ======================================================================

template<int W>
class FlatHash {
  int capacity_;
  int count_;
  std::vector<uint64_t> data_;
  std::vector<uint8_t> status_;

  size_t slot_for(const uint64_t* bs) const {
    return BitOps<W>::hash(bs) & static_cast<size_t>(capacity_ - 1);
  }

  bool insert_internal(const uint64_t* bs) {
    size_t slot = slot_for(bs);
    while (true) {
      if (status_[slot] == 0) {
        BitOps<W>::copy(&data_[slot * W], bs);
        status_[slot] = 1;
        count_++;
        return true;
      }
      if (BitOps<W>::equal(&data_[slot * W], bs)) return false;
      slot = (slot + 1) & static_cast<size_t>(capacity_ - 1);
    }
  }

  void grow() {
    int old_cap = capacity_;
    std::vector<uint64_t> old_data;
    std::vector<uint8_t> old_status;
    old_data.swap(data_);
    old_status.swap(status_);

    capacity_ *= 2;
    data_.assign(static_cast<size_t>(capacity_) * W, 0);
    status_.assign(capacity_, 0);
    count_ = 0;

    for (int i = 0; i < old_cap; i++) {
      if (old_status[i] == 1) {
        insert_internal(&old_data[static_cast<size_t>(i) * W]);
      }
    }
  }

public:
  explicit FlatHash(int initial_capacity = 1024) : count_(0) {
    capacity_ = 1;
    while (capacity_ < initial_capacity) capacity_ *= 2;
    data_.assign(static_cast<size_t>(capacity_) * W, 0);
    status_.assign(capacity_, 0);
  }

  bool insert(const uint64_t* bs) {
    if (count_ * 10 >= capacity_ * 7) grow();
    return insert_internal(bs);
  }

  bool contains(const uint64_t* bs) const {
    size_t slot = slot_for(bs);
    while (true) {
      if (status_[slot] == 0) return false;
      if (BitOps<W>::equal(&data_[slot * W], bs)) return true;
      slot = (slot + 1) & static_cast<size_t>(capacity_ - 1);
    }
  }

  int size() const { return count_; }
};


// ======================================================================
// Section 4: Decode helper (bitset -> R integer vector) and full set
// ======================================================================

inline Rcpp::IntegerVector decode_bitset(const uint64_t* bs, int W, int n) {
  std::vector<int> result;
  for (int i = 0; i < n; i++) {
    if (bs_test_bit(bs, i)) result.push_back(i + 1);
  }
  return Rcpp::wrap(result);
}

inline void compute_full_set(uint64_t* fs, int W, int n) {
  for (int i = 0; i < W - 1; i++) fs[i] = ~0ULL;
  int remainder = n % 64;
  fs[W - 1] = (remainder == 0) ? ~0ULL : ((1ULL << remainder) - 1);
}


// ======================================================================
// Section 5: Runtime-W bitset operations (for W > 3)
// ======================================================================

namespace rtops {
  inline void bor(uint64_t* d, const uint64_t* a, const uint64_t* b, int W) {
    for (int i = 0; i < W; i++) d[i] = a[i] | b[i]; }
  inline void band(uint64_t* d, const uint64_t* a, const uint64_t* b, int W) {
    for (int i = 0; i < W; i++) d[i] = a[i] & b[i]; }
  inline void bnot(uint64_t* d, const uint64_t* a, const uint64_t* m, int W) {
    for (int i = 0; i < W; i++) d[i] = ~a[i] & m[i]; }
  inline bool equal(const uint64_t* a, const uint64_t* b, int W) {
    return std::memcmp(a, b, W * sizeof(uint64_t)) == 0; }
  inline bool is_zero(const uint64_t* a, int W) {
    for (int i = 0; i < W; i++) {
      if (a[i]) return false;
    }
    return true;
  }
  inline bool is_subset(const uint64_t* a, const uint64_t* b, int W) {
    for (int i = 0; i < W; i++) {
      if ((a[i] & b[i]) != a[i]) return false;
    }
    return true;
  }
  inline size_t hash(const uint64_t* a, int W) {
    size_t h = 0;
    for (int i = 0; i < W; i++) {
      uint64_t x = a[i]; x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
      x ^= x >> 27; x *= 0x94d049bb133111ebULL; x ^= x >> 31;
      h ^= static_cast<size_t>(x) + 0x9e3779b9 + (h << 6) + (h >> 2);
    }
    return h;
  }
}


// ======================================================================
// Section 6: Runtime-W FlatVec and FlatHash
// ======================================================================

struct RtVec {
  int W; std::vector<uint64_t> data;
  explicit RtVec(int w) : W(w) {}
  int size() const { return static_cast<int>(data.size()) / W; }
  uint64_t* get(int i) { return &data[static_cast<size_t>(i) * W]; }
  const uint64_t* get(int i) const { return &data[static_cast<size_t>(i) * W]; }
  void push(const uint64_t* bs) { data.insert(data.end(), bs, bs + W); }
  void reserve(int c) { data.reserve(static_cast<size_t>(c) * W); }
};

class RtHash {
  int W_, capacity_, count_;
  std::vector<uint64_t> data_; std::vector<uint8_t> status_;
  size_t slot_for(const uint64_t* bs) const {
    return rtops::hash(bs, W_) & static_cast<size_t>(capacity_ - 1); }
  bool insert_internal(const uint64_t* bs) {
    size_t s = slot_for(bs);
    while (true) {
      if (!status_[s]) { std::memcpy(&data_[s*W_], bs, W_*8); status_[s]=1; count_++; return true; }
      if (rtops::equal(&data_[s*W_], bs, W_)) return false;
      s = (s+1) & static_cast<size_t>(capacity_-1);
    }
  }
  void grow() {
    int oc = capacity_; std::vector<uint64_t> od; std::vector<uint8_t> os;
    od.swap(data_); os.swap(status_);
    capacity_ *= 2; data_.assign(static_cast<size_t>(capacity_)*W_, 0);
    status_.assign(capacity_, 0); count_ = 0;
    for (int i = 0; i < oc; i++) if (os[i]) insert_internal(&od[static_cast<size_t>(i)*W_]);
  }
public:
  RtHash(int w, int ic = 1024) : W_(w), count_(0) {
    capacity_ = 1; while (capacity_ < ic) capacity_ *= 2;
    data_.assign(static_cast<size_t>(capacity_)*W_, 0); status_.assign(capacity_, 0);
  }
  bool insert(const uint64_t* bs) { if (count_*10 >= capacity_*7) grow(); return insert_internal(bs); }
  bool contains(const uint64_t* bs) const {
    size_t s = slot_for(bs);
    while (true) { if (!status_[s]) return false; if (rtops::equal(&data_[s*W_], bs, W_)) return true;
      s = (s+1) & static_cast<size_t>(capacity_-1); }
  }
  int size() const { return count_; }
};


#endif // TOPOLOGYR_BITSET_OPS_H_
