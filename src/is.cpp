#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <climits>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <utility>
#include <vector>

using namespace std;

typedef struct {
  unsigned height;
  int width;
  unsigned client;
} item;

typedef struct {
  std::pair<int, int> bottom_point;
  int width;
} mos_t;

struct bottom_left_cmp_open_space {
  bool operator()(const mos_t &a, const mos_t &b) const
  {
    if (a.bottom_point.second != b.bottom_point.second) {
      return a.bottom_point.second < b.bottom_point.second;
    }
    return a.bottom_point.first < b.bottom_point.first;
  }
};

struct dominated_cmp {
  bool operator()(const mos_t &a, const mos_t &b) const
  {
    if (a.bottom_point.first != b.bottom_point.first) {
      return a.bottom_point.first < b.bottom_point.first;
    }
    return a.bottom_point.second < b.bottom_point.second;
  }
};

struct SegmentTree {
  struct Node {
    int l, r;
    int min_val;
    bool lazy;
    int lazy_val;
    Node *left, *right;

    Node(int l, int r)
        : l(l),
          r(r),
          min_val(INT_MAX),
          lazy(false),
          lazy_val(INT_MAX),
          left(nullptr),
          right(nullptr)
    {
    }
  };

  Node *root;

  SegmentTree(int n) { root = build(0, n - 1); }

  ~SegmentTree() { destroy(root); }

  void destroy(Node *node)
  {
    if (!node) return;
    destroy(node->left);
    destroy(node->right);
    delete node;
  }

  Node *build(int l, int r)
  {
    Node *node = new Node(l, r);
    if (l < r) {
      int mid = (l + r) / 2;
      node->left = build(l, mid);
      node->right = build(mid + 1, r);
    }
    return node;
  }

  void apply(Node *node, int val)
  {
    node->min_val = val;
    node->lazy = true;
    node->lazy_val = val;
  }

  void push(Node *node)
  {
    if (node->lazy && node->left) {
      apply(node->left, node->lazy_val);
      apply(node->right, node->lazy_val);
      node->lazy = false;
    }
  }

  void update(Node *node, int ql, int qr, int val)
  {
    if (qr < node->l || node->r < ql) return;

    if (ql <= node->l && node->r <= qr) {
      apply(node, val);
      return;
    }

    if (node->lazy) {
      push(node);
    }

    int mid = (node->l + node->r) / 2;
    if (ql <= mid) update(node->left, ql, qr, val);
    if (mid + 1 <= qr) update(node->right, ql, qr, val);

    node->min_val = std::min(node->left->min_val, node->right->min_val);
    node->lazy = false;
  }

  int query(Node *node, int ql, int qr)
  {
    if (qr < node->l || node->r < ql) return INT_MAX;

    if (ql <= node->l && node->r <= qr) {
      return node->min_val;
    }

    if (node->lazy) {
      return node->min_val;
    }

    int mid = (node->l + node->r) / 2;
    int leftMin = INT_MAX, rightMin = INT_MAX;

    if (ql <= mid) leftMin = query(node->left, ql, qr);
    if (mid + 1 <= qr) rightMin = query(node->right, ql, qr);

    return std::min(leftMin, rightMin);
  }

  void update(int l, int r, int val) { update(root, l, r, val); }
  int query(int l, int r) { return query(root, l, r); }
};

inline bool can_item_fit(const item &item, const mos_t &mos)
{
  return item.width <= mos.width;
}

void update_strip_height(unsigned &strip_height, int &mos_index,
                         const unsigned &item_height,
                         const std::vector<mos_t> &slist)
{
  mos_t mos = slist[mos_index];
  unsigned height = mos.bottom_point.second + item_height;
  if (height > strip_height) {
    strip_height = height;
  }
}

bool intersects(const mos_t &mos, int space_bp_x, int space_bp_y,
                unsigned space_width, const unsigned &item_width,
                const unsigned &item_height)
{
  int item_x1 = space_bp_x;
  int item_x2 = space_bp_x + (int)item_width;
  int item_y1 = space_bp_y;
  int item_y2 = space_bp_y + (int)item_height;

  int mos_x1 = mos.bottom_point.first;
  int mos_x2 = mos_x1 + mos.width;
  int mos_y1 = mos.bottom_point.second;
  int mos_y2 = std::numeric_limits<int>::max();

  return (item_x1 < mos_x2 && item_x2 > mos_x1 && item_y1 < mos_y2 &&
          item_y2 > mos_y1);
}

bool is_dominated(const mos_t &p, const std::vector<mos_t> &layer)
{
  auto range = std::equal_range(
      layer.begin(), layer.end(), p, [](const mos_t &a, const mos_t &b) {
        return a.bottom_point.first < b.bottom_point.first;
      });

  for (auto it = range.first; it != range.second; ++it) {
    if (it->bottom_point.second <= p.bottom_point.second &&
        it->width >= p.width) {
      return true;
    }
  }
  return false;
}

void place_item(const item &item, int &mos_index, std::vector<mos_t> &slist,
                bool debug_sol = false, fstream *solfile = nullptr)
{
  mos_t space = slist[mos_index];
  int space_width = space.width;
  int item_tp_x = space.bottom_point.first + item.width;
  int space_bp_x = space.bottom_point.first;
  int space_bp_y = space.bottom_point.second;

  if (debug_sol) {
    *solfile << space_bp_x << " " << space.bottom_point.second << "\n"
             << item_tp_x << " " << space.bottom_point.second + item.height
             << "\n";
  }

  slist.erase(slist.begin() + mos_index);
  std::vector<mos_t> to_process;

  mos_t new_space1 = {{space_bp_x, space_bp_y + item.height}, space_width};
  to_process.push_back(new_space1);

  if (space_width - item.width > 0) {
    mos_t new_space2 = {{space_bp_x + item.width, space_bp_y},
                        space_width - item.width};
    to_process.push_back(new_space2);
  }

  for (auto mos = slist.begin(); mos != slist.end();) {
    if (intersects(*mos, space_bp_x, space_bp_y, space_width, item.width,
                   item.height)) {
      int mos_x = mos->bottom_point.first;
      int mos_y = mos->bottom_point.second;
      int mos_width = mos->width;

      mos = slist.erase(mos);

      if (mos_x < space_bp_x) {
        mos_t new_space3 = {{mos_x, mos_y}, space_bp_x - mos_x};
        to_process.push_back(new_space3);
      }

      if (mos_x + mos_width > item_tp_x) {
        mos_t new_space4 = {{item_tp_x, mos_y}, mos_x + mos_width - item_tp_x};
        to_process.push_back(new_space4);
      }

      mos_t new_space5 = {{mos_x, space_bp_y + item.height}, mos_width};
      to_process.push_back(new_space5);
    }
    else {
      mos++;
    }
  }

  std::vector<mos_t> slistByX = slist;
  std::sort(slistByX.begin(), slistByX.end(), dominated_cmp());

  for (auto new_space : to_process) {
    if (!is_dominated(new_space, slistByX)) {
      slist.push_back(new_space);
      slistByX.push_back(new_space);
      std::sort(slistByX.begin(), slistByX.end(), dominated_cmp());
    }
  }

  std::sort(slist.begin(), slist.end(), bottom_left_cmp_open_space());
}

int find_mos(const item &item, std::vector<mos_t> &slist, SegmentTree &seg_tree)
{
  for (size_t i = 0; i < slist.size(); i++) {
    if (can_item_fit(item, slist[i])) {
      int space_bp_x = slist[i].bottom_point.first;
      int space_bp_y = slist[i].bottom_point.second;
      unsigned item_tp_x = slist[i].bottom_point.first + item.width;
      int min_client = seg_tree.query(space_bp_x, item_tp_x - 1);
      if (min_client >= item.client) {
        seg_tree.update(space_bp_x, item_tp_x - 1, item.client);
        return i;
      }
    }
  }
  return -1;
}

unsigned first_fit_pack(const int &max_width, const std::vector<item> &items,
                        const std::vector<unsigned> &sequence,
                        bool debug_sol = false, fstream *solfile = nullptr)
{
  unsigned strip_height = 0;
  std::vector<mos_t> slist = {{{0, 0}, max_width}};
  SegmentTree seg_tree(max_width);

  unsigned items_count = sequence.size();

  for (unsigned i = 0; i < items_count; i++) {
    unsigned item_index = sequence[i];
    const item &item = items[item_index];

    int mos_index = find_mos(item, slist, seg_tree);
    if (mos_index == -1) {
      return std::numeric_limits<unsigned>::max();
    }

    update_strip_height(strip_height, mos_index, item.height, slist);
    place_item(item, mos_index, slist, debug_sol, solfile);

    if (debug_sol) {
      *solfile << item_index << "\n" << item.client << "\n";
    }
  }

  return strip_height;
}

void swap_random(std::vector<unsigned> &seq_copy, std::mt19937 &rng)
{
  std::uniform_int_distribution<size_t> dist(0, seq_copy.size() - 1);
  size_t i = dist(rng);
  size_t j = dist(rng);
  while (j == i) j = dist(rng);
  std::swap(seq_copy[i], seq_copy[j]);
}

std::pair<unsigned, std::vector<unsigned>> random_ls(
    const unsigned &max_width, const std::vector<item> &items,
    const std::vector<unsigned> &sequence, unsigned iter, std::mt19937 &rng,
    const std::chrono::high_resolution_clock::time_point &start,
    int time_limit_s, bool &timeout)
{
  std::vector<unsigned> best_seq = sequence;
  unsigned height = first_fit_pack(max_width, items, sequence);

  using clock = std::chrono::high_resolution_clock;

  for (unsigned i = 1; i <= iter; i++) {
    auto now = clock::now();
    auto elapsed =
        std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
    if (elapsed >= time_limit_s) {
      timeout = true;
      break;
    }

    std::vector<unsigned> new_seq = best_seq;
    swap_random(new_seq, rng);
    unsigned new_height = first_fit_pack(max_width, items, new_seq);

    if (new_height <= height) {
      height = new_height;
      best_seq = new_seq;
    }
  }

  return {height, best_seq};
}

bool sort_by_client_and_area_desc(unsigned a, unsigned b,
                                  const std::vector<item> &items)
{
  if (items[a].client != items[b].client) {
    return items[a].client > items[b].client;
  }
  return (items[a].height * items[a].width) >
         (items[b].height * items[b].width);
}

bool sort_by_client_and_width_desc(unsigned a, unsigned b,
                                   const std::vector<item> &items)
{
  if (items[a].client != items[b].client) {
    return items[a].client > items[b].client;
  }
  return items[a].width > items[b].width;
}

bool sort_by_client_and_height_desc(unsigned a, unsigned b,
                                    const std::vector<item> &items)
{
  if (items[a].client != items[b].client) {
    return items[a].client > items[b].client;
  }
  return items[a].height > items[b].height;
}

void build_initial_seqs(std::vector<unsigned> seq,
                        const std::vector<item> &items,
                        std::vector<std::vector<unsigned>> &initial_seqs)
{
  initial_seqs.resize(3);

  std::sort(seq.begin(), seq.end(), [&items](unsigned a, unsigned b) {
    return sort_by_client_and_area_desc(a, b, items);
  });
  initial_seqs[0] = seq;

  std::sort(seq.begin(), seq.end(), [&items](unsigned a, unsigned b) {
    return sort_by_client_and_width_desc(a, b, items);
  });
  initial_seqs[1] = seq;

  std::sort(seq.begin(), seq.end(), [&items](unsigned a, unsigned b) {
    return sort_by_client_and_height_desc(a, b, items);
  });
  initial_seqs[2] = seq;
}

unsigned iterative_search(const int &time_limit_s, const unsigned &max_width,
                          const std::vector<item> &items,
                          std::vector<unsigned> &sequence, std::mt19937 &rng)
{
  using clock = std::chrono::high_resolution_clock;
  auto start = clock::now();
  unsigned min_height = std::numeric_limits<unsigned>::max();
  unsigned iter = 1;

  std::vector<std::vector<unsigned>> initial_seqs;
  build_initial_seqs(sequence, items, initial_seqs);

  bool timeout = false;

  while (!timeout) {
    for (auto &seq : initial_seqs) {
      auto [height, best_sequence] = random_ls(max_width, items, seq, iter, rng,
                                               start, time_limit_s, timeout);

      if (height < min_height) {
        min_height = height;
        sequence = best_sequence;
      }

      if (timeout) break;
    }
    iter *= 2;
  }

  return min_height;
}

void show_help(const char *progname)
{
  cout << "Usage: " << progname << " [options]\n"
       << "Options:\n"
       << "  -f <file>      Specifies the path of the input instance\n"
       << "  -s <string>    Sets the seed for generating random numbers\n"
       << "  -t <number>    Sets the time limit for execution (defaults to "
          "60s)\n"
       << "  -d             Activates debug mode, the output is a file you can "
          "use as input in debug.py script\n"
       << "  -h             Shows this help menu\n";
}

int main(int argc, char *argv[])
{
  int opt;
  string input_filename;
  unsigned seed = 1;
  int time_limit_s = 60;
  bool debug_sol = false;

  while ((opt = getopt(argc, argv, "f:s:t:dh")) != -1) {
    switch (opt) {
      case 'f':
        input_filename = optarg;
        break;

      case 's':
        seed = atoi(optarg);
        break;

      case 't':
        time_limit_s = atoi(optarg);
        break;

      case 'd':
        debug_sol = true;
        break;

      case 'h':
        show_help(argv[0]);
        return 0;

      default:
        show_help(argv[0]);
        return 1;
    }
  }

  ifstream input_file(input_filename);

  unsigned num_items;
  int max_width;
  unsigned num_clients, num_items_client, client, width, height;

  input_file >> num_clients >> num_items >> max_width;

  vector<item> items(num_items);
  vector<unsigned> seq(num_items);

  int index = 0;

  for (unsigned i = 0; i < num_clients; i++) {
    input_file >> client;
    input_file >> num_items_client;

    for (unsigned j = 0; j < num_items_client; j++) {
      input_file >> height >> width;
      items[index].height = height;
      items[index].width = width;
      items[index].client = client;

      seq[index] = index;

      index++;
    }
  }

  int val[num_items] = {0};
  int area = 0;
  int lb2 = 0;
  int lb1, lb;
  unsigned ub = 0;

  for (unsigned i = 0; i < num_items; i++) {
    ub += items[i].height;
    area += items[i].width * items[i].height;
    val[i] += items[i].height;
    lb2 = max(lb2, val[i]);

    for (unsigned j = i + 1; j < num_items; j++) {
      if (items[j].client != items[i].client &&
          items[i].width + items[j].width > max_width && val[i] > val[j]) {
        val[j] = val[i];
      }
    }
  }

  lb1 = area / max_width;
  lb = max(lb1, lb2);

  const std::string logs_folder = "logs";

  if (!std::filesystem::exists(logs_folder)) {
    std::filesystem::create_directory(logs_folder);
  }

  std::filesystem::path filePath(input_filename);
  std::filesystem::path solPath(input_filename);
  std::string parentPathStr = filePath.parent_path().string();
  std::size_t pos = parentPathStr.find("instances/");
  std::string relativePath =
      parentPathStr.substr(pos + std::string("instances/").length());
  std::filesystem::path logFilePath = std::filesystem::path(logs_folder) /
                                      relativePath /
                                      (filePath.stem().string() + ".log");

  if (!std::filesystem::exists(logFilePath.parent_path())) {
    std::filesystem::create_directories(logFilePath.parent_path());
  }

  std::fstream logfile(logFilePath, std::ios::app);

  std::mt19937 rng(seed);

  using clock = std::chrono::high_resolution_clock;
  auto start = clock::now();

  unsigned min_height =
      iterative_search(time_limit_s, max_width, items, seq, rng);

  auto end = clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start)
          .count();

  if (debug_sol) {
    const std::string sol_folder = "sol";

    if (!std::filesystem::exists(sol_folder)) {
      std::filesystem::create_directory(sol_folder);
    }

    std::filesystem::path solFilePath = std::filesystem::path(sol_folder) /
                                        relativePath /
                                        (filePath.stem().string() + "Sol.txt");

    if (!std::filesystem::exists(solFilePath.parent_path())) {
      std::filesystem::create_directories(solFilePath.parent_path());
    }

    std::fstream solfile(solFilePath, std::ios::out);

    solfile << max_width << "\n";
    solfile << min_height << "\n";
    int dummy = first_fit_pack(max_width, items, seq, debug_sol, &solfile);
  }

  cout << "Minimum height: " << min_height << endl;
  cout << (double)min_height / lb << endl;

  cout << "Tempo de execução: " << (double)duration / 1000000.0 << "s\n";

  logfile << "Largura da faixa = " << max_width << "\n";
  logfile << "Melhor altura = " << min_height << "\n\n";
  logfile << "Lower Bound 1 = " << lb1 << "\n";
  logfile << "Lower Bound 2 = " << lb2 << "\n";
  logfile << "Avaliação = " << (double)min_height / lb << "\n\n";

  logfile << "Tempo de execução: " << (double)duration / 1000000.0 << "s"
          << "\n\n"
          << "---------------------------" << "\n\n";

  return 0;
}