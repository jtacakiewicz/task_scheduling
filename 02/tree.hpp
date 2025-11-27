#include <vector>


template<class T>
class IntervalTree {
    struct Interval {
        T id;
        long long l, r; // [l, r]
    };
    struct Node {
        Interval iv;
        long long max_r; // maximum r in this subtree
        Node *left, *right;
        Node(const Interval& iv): iv(iv), max_r(iv.r), left(nullptr), right(nullptr) {}
    };
    Node* root = nullptr;

public:
    IntervalTree() = default;

    void insert(T id, long long l, long long r) {
        root = insert(root, Interval{id, l, r});
    }

    // return all intervals that overlap [ql, qr]
    std::vector<Interval*> query(long long ql, long long qr) const {
        std::vector<Interval*> res;
        query(root, ql, qr, res);
        return res;
    }
private:
    static bool overlap(const Interval& a, long long ql, long long qr) {
        return !(a.r < ql || qr < a.l);
    }
    Node* insert(Node* n, Interval iv) {
        if (!n) return new Node(iv);
        if (iv.l < n->iv.l)
            n->left = insert(n->left, iv);
        else
            n->right = insert(n->right, iv);
        n->max_r = n->iv.r;
        if (n->left)  n->max_r = std::max(n->max_r, n->left->max_r);
        if (n->right) n->max_r = std::max(n->max_r, n->right->max_r);
        return n;
    }
    void query(Node* n, long long ql, long long qr, std::vector<Interval*>& out) const {
        if (!n) return;
        if (n->max_r < ql) return;
        // left subtree first
        query(n->left, ql, qr, out);
        // check this interval
        if (overlap(n->iv, ql, qr)) out.push_back(&n->iv);
        // explore right subtree only if it may overlap
        if (n->iv.l <= qr)
            query(n->right, ql, qr, out);
    }
};
