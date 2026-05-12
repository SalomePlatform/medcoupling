// Copyright (C) 2007-2026  CEA, EDF
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#pragma once

#include "InterpKernelException.hxx"

#include <vector>
#include <algorithm>

#include <iostream>
#include <stack>
#include <memory>
#include <limits>
#include <cmath>

constexpr double BBTREE_DFT_EPSILON = 1e-12;

template <int dim, class ConnType, class Derived>
class BBTreeBaseBaseBase
{
   protected:
    std::unique_ptr<Derived> _left;

    std::unique_ptr<Derived> _right;

    int _level = 0;

    const double *_bb = nullptr;

   protected:
    bool terminal() const { return !_left && !_right; }

    Derived *left() const { return _left.get(); }

    Derived *right() const { return _right.get(); }

   public:
    const double *getBBoxData() const { return this->_bb; }

    template <class TerminalVertexOwner>
    void walk(std::vector<TerminalVertexOwner> &container)
    {
        if (this->terminal())
        {
            if (!static_cast<Derived *>(this)->empty())
            {
                container.push_back(TerminalVertexOwner(static_cast<Derived *>(this)));
            }
            return;
        }
        this->_left->walk(container);
        this->_right->walk(container);
    }

   protected:
    BBTreeBaseBaseBase() = default;

    BBTreeBaseBaseBase(int level) : _level(level) {}

    BBTreeBaseBaseBase(const double *bbs, int level) : _level(level), _bb(bbs) {}

    static std::unique_ptr<double[]> FillDataNotTerminal(
        int level, const double *bbs, const ConnType *elems, ConnType nbelems
    )
    {
        std::unique_ptr<double[]> nodes(new double[nbelems]);
        for (ConnType i = 0; i < nbelems; i++)
        {
            ConnType elem;
            if (elems)
                elem = elems[i];
            else
                elem = i;
            nodes[i] = bbs[elem * dim * 2 + (level % dim) * 2];
        }
        return nodes;
    }

    static void SplitForNotTerminal(
        int level,
        const double *bbs,
        const ConnType *elems,
        ConnType nbelems,
        std::vector<ConnType> &left,
        double &max_left,
        std::vector<ConnType> &right,
        double &min_right
    )
    {
        std::unique_ptr<double[]> nodes(FillDataNotTerminal(level, bbs, elems, nbelems));
        std::nth_element<double *>(nodes.get(), nodes.get() + nbelems / 2, nodes.get() + nbelems);
        double median(*(nodes.get() + nbelems / 2));

        left.reserve(nbelems / 2 + 1);
        right.reserve(nbelems / 2 + 1);
        max_left = -std::numeric_limits<double>::max();
        min_right = std::numeric_limits<double>::max();
        for (ConnType i = 0; i < nbelems; i++)
        {
            ConnType elem(elems ? elems[i] : i);
            double max(bbs[elem * dim * 2 + (level % dim) * 2 + 1]);
            double min(bbs[elem * dim * 2 + (level % dim) * 2]);

            if (min > median)
            {
                right.push_back(elem);
                if (min < min_right)
                    min_right = min;
            }
            else
            {
                left.push_back(elem);
                if (max > max_left)
                    max_left = max;
            }
        }
    }

   public:
    // DFS iterator implementation to allow to iterate across terminal leaves
    class LeafIterator
    {
       private:
        std::stack<const BBTreeBaseBaseBase<dim, ConnType, Derived> *> stack;
        const BBTreeBaseBaseBase<dim, ConnType, Derived> *_current = nullptr;

        void advanceToNextLeaf()
        {
            _current = nullptr;

            while (!stack.empty())
            {
                const BBTreeBaseBaseBase<dim, ConnType, Derived> *node(stack.top());
                stack.pop();
                if (!node)
                {
                    continue;
                }

                if (node->terminal())
                {
                    _current = node;
                    return;
                }
                // DFS
                if (node->right())
                    stack.push(node->right());
                if (node->left())
                    stack.push(node->left());
            }
        }

       public:
        using iterator_category = std::input_iterator_tag;
        using value_type = BBTreeBaseBaseBase<dim, ConnType, Derived>;
        using difference_type = std::ptrdiff_t;
        using pointer = const BBTreeBaseBaseBase<dim, ConnType, Derived> *;
        using reference = const BBTreeBaseBaseBase<dim, ConnType, Derived> &;

        LeafIterator() = default;

        explicit LeafIterator(const BBTreeBaseBaseBase<dim, ConnType, Derived> *root)
        {
            if (root)
            {
                stack.push(root);
            }
            advanceToNextLeaf();
        }

        reference operator*() const { return *_current; }
        pointer operator->() const { return _current; }

        LeafIterator &operator++()
        {
            advanceToNextLeaf();
            return *this;
        }

        bool operator==(const LeafIterator &other) const { return _current == other._current; }

        bool operator!=(const LeafIterator &other) const { return !(*this == other); }
    };

    LeafIterator begin() const { return LeafIterator(this); }
    LeafIterator end() const { return LeafIterator(); }
    // End of DFS iterator implementation to allow to iterate across terminal leaves
};

template <int dim, class ConnType, int MIN_NB_ELEMS, int MAX_LEVEL, class Derived>
class BBTreeBaseBase : public BBTreeBaseBaseBase<dim, ConnType, Derived>
{
   protected:
    typename std::vector<ConnType> _elems;

   protected:
    constexpr bool constexprTerminal(int level, ConnType nbelems)
    {
        return nbelems < MIN_NB_ELEMS || level > MAX_LEVEL;
    }

    ConnType size()
    {
        if (this->terminal())
            return (ConnType)this->_elems.size();
        return this->_left->size() + this->_right->size();
    }

   public:
    const std::vector<ConnType> &getElements() const { return this->_elems; }

    static void FillDataTerminal(const ConnType *elems, ConnType nbelems, typename std::vector<ConnType> &elemsOut)
    {
        elemsOut.resize(nbelems);
        for (ConnType i = 0; i < nbelems; i++)
        {
            ConnType elem;
            if (elems)
                elem = elems[i];
            else
                elem = i;

            elemsOut[i] = elem;
        }
    }

    // required for walk concept
    bool empty() const { return this->_elems.empty(); }

   protected:
    BBTreeBaseBase() = default;

    BBTreeBaseBase(int level) : BBTreeBaseBaseBase<dim, ConnType, Derived>(level) {}

    BBTreeBaseBase(const double *bbs, const ConnType *elems, int level, ConnType nbelems)
        : BBTreeBaseBaseBase<dim, ConnType, Derived>(bbs, level)
    {
        if (constexprTerminal(level, nbelems))
        {
            FillDataTerminal(elems, nbelems, _elems);
            return;
        }
    }
};

template <int dim, class ConnType, int MIN_NB_ELEMS, int MAX_LEVEL, class Derived>
class BBTreeBase : public BBTreeBaseBase<dim, ConnType, MIN_NB_ELEMS, MAX_LEVEL, Derived>
{
   private:
    double _max_left;
    double _min_right;
    double _epsilon;

   public:
    /*!
      Constructor of the bounding box tree
      \param bbs pointer to the [xmin1 xmax1 ymin1 ymax1 xmin2 xmax2 ...] array containing the bounding boxes that are
      to be indexed.
      \param elems array to the indices of the elements contained in the BBTree
      \param level level in the BBTree recursive structure
      \param nbelems nb of elements in the BBTree
      \param epsilon precision to which points are decided to be coincident. Epsilon can be positive or negative.
             If \a epsilon is positive the request method will enlarge the computed bounding box (more matching elems
      return). If negative the given bounding box will be tighten (less matching elems return).

      Parameters \a elems and \a level are used only by BBTree itself for creating trees recursively. A typical use is
      therefore :
      \code
      int nbelems=...
      double* bbs= new double[2*2*nbelems];
      // filling bbs ...
      ...
      BBTree<2> tree = new BBTree<2>(elems,0,0,nbelems,1e-12);
      \endcode
    */
    BBTreeBase(const double *bbs, const ConnType *elems, int level, ConnType nbelems, double epsilon)
        : BBTreeBaseBase<dim, ConnType, MIN_NB_ELEMS, MAX_LEVEL, Derived>(bbs, elems, level, nbelems), _epsilon(epsilon)
    {
        if (this->constexprTerminal(level, nbelems))
        {
            return;
        }
        std::vector<ConnType> new_elems_left, new_elems_right;
        double max_left, min_right;
        this->SplitForNotTerminal(level, bbs, elems, nbelems, new_elems_left, max_left, new_elems_right, min_right);
        _max_left = max_left + std::abs(_epsilon);
        _min_right = min_right - std::abs(_epsilon);
        auto FromVectToPt = [](const std::vector<ConnType> &v) -> const ConnType *
        { return v.empty() ? nullptr : v.data(); };
        this->_left.reset(
            new Derived(bbs, FromVectToPt(new_elems_left), level + 1, (ConnType)new_elems_left.size(), _epsilon)
        );
        this->_right.reset(
            new Derived(bbs, FromVectToPt(new_elems_right), level + 1, (ConnType)new_elems_right.size(), _epsilon)
        );
    }

    /*! returns in \a elems the list of elements potentially intersecting the bounding box pointed to by \a bb

      \param bb pointer to query bounding box
      \param elems list of elements (given in 0-indexing that is to say in \b C \b mode) intersecting the bounding box
    */
    void getIntersectingElems(const double *bb, std::vector<ConnType> &elems) const
    {
        //  terminal node : return list of elements intersecting bb
        if (this->terminal())
        {
            for (std::size_t i = 0; i < this->_elems.size(); i++)
            {
                const double *const bb_ptr = this->_bb + this->_elems[i] * 2 * dim;
                bool intersects = true;
                for (int idim = 0; idim < dim; idim++)
                {
                    if (bb_ptr[idim * 2] - bb[idim * 2 + 1] > -_epsilon ||
                        bb_ptr[idim * 2 + 1] - bb[idim * 2] < _epsilon)
                        intersects = false;
                }
                if (intersects)
                {
                    elems.push_back(this->_elems[i]);
                }
            }
            return;
        }

        // non terminal node
        double min = bb[(this->_level % dim) * 2];
        double max = bb[(this->_level % dim) * 2 + 1];
        if (max < _min_right)
        {
            this->_left->getIntersectingElems(bb, elems);
            return;
        }
        if (min > _max_left)
        {
            this->_right->getIntersectingElems(bb, elems);
            return;
        }
        this->_left->getIntersectingElems(bb, elems);
        this->_right->getIntersectingElems(bb, elems);
    }

    /*!
     * This method is very close to getIntersectingElems except that it returns number of elems instead of elems
     * themselves.
     */
    ConnType getNbOfIntersectingElems(const double *bb)
    {
        //  terminal node : return list of elements intersecting bb
        ConnType ret(0);
        if (this->terminal())
        {
            for (std::size_t i = 0; i < this->_elems.size(); i++)
            {
                const double *const bb_ptr = this->_bb + this->_elems[i] * 2 * dim;
                bool intersects = true;
                for (int idim = 0; idim < dim; idim++)
                {
                    if (bb_ptr[idim * 2] - bb[idim * 2 + 1] > -_epsilon ||
                        bb_ptr[idim * 2 + 1] - bb[idim * 2] < _epsilon)
                        intersects = false;
                }
                if (intersects)
                    ret++;
            }
            return ret;
        }
        // non terminal node
        double min = bb[(this->_level % dim) * 2];
        double max = bb[(this->_level % dim) * 2 + 1];
        if (max < _min_right)
            return this->_left->getNbOfIntersectingElems(bb);
        if (min > _max_left)
            return this->_right->getNbOfIntersectingElems(bb);
        return this->_left->getNbOfIntersectingElems(bb) + this->_right->getNbOfIntersectingElems(bb);
    }

    /*! returns in \a elems the list of elements potentially containing the point pointed to by \a xx
      \param xx pointer to query point coords
      \param elems list of elements (given in 0-indexing) intersecting the bounding box
    */
    void getElementsAroundPoint(const double *xx, std::vector<ConnType> &elems) const
    {
        //  terminal node : return list of elements intersecting bb
        if (this->terminal())
        {
            for (std::size_t i = 0; i < this->_elems.size(); i++)
            {
                const double *const bb_ptr = this->_bb + this->_elems[i] * 2 * dim;
                bool intersects = true;
                for (int idim = 0; idim < dim; idim++)
                {
                    if (bb_ptr[idim * 2] - xx[idim] > _epsilon || bb_ptr[idim * 2 + 1] - xx[idim] < -_epsilon)
                        intersects = false;
                }
                if (intersects)
                {
                    elems.push_back(this->_elems[i]);
                }
            }
            return;
        }

        // non terminal node
        if (xx[this->_level % dim] < _min_right)
        {
            this->_left->getElementsAroundPoint(xx, elems);
            return;
        }
        if (xx[this->_level % dim] > _max_left)
        {
            this->_right->getElementsAroundPoint(xx, elems);
            return;
        }
        this->_left->getElementsAroundPoint(xx, elems);
        this->_right->getElementsAroundPoint(xx, elems);
    }
};

template <int dim, class ConnType = int>
class BBTree : public BBTreeBase<dim, ConnType, 15, 20, BBTree<dim, ConnType> >
{
   public:
    BBTree(const double *bbs, const ConnType *elems, int level, ConnType nbelems, double epsilon = BBTREE_DFT_EPSILON)
        : BBTreeBase<dim, ConnType, 15, 20, BBTree<dim, ConnType> >(bbs, elems, level, nbelems, epsilon)
    {
    }
};
