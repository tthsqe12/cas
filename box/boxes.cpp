#include "boxes.h"
#include "ex_parse.h"
#include <stack>
#include "serialize.h"
#include "notebook.h"
#include "box_convert.h"

/*
const char* cell_type_names[cellt_MAX] = {
    "Input",
    "Output",
    "Message",
    "Print",
    "MathCaption",
    "Caption",
    "ExampleText",
    "Text",
    "Definition",
    "Subsubsection",
    "Subsection",
    "Section",
    "Title"
};
*/


std::string stdvector_tostring(std::vector<int32_t> v)
{
    std::string s = "<";
    for (size_t i = 0; i < v.size(); i++)
    {
        s.append(stdstring_to_string(v[i]));
        if (i + 1 < v.size())
            s.push_back(',');
    }
    s.push_back('>');
    return s;
}


ibox iboxchar_set_sizes(ibox b, double sizex_, double sizey_, double centery_)
{
    int32_t sizex = 0.5 + IMM_UNITS_PER_PIXEL*sizex_;
    int32_t sizey = 0.5 + IMM_UNITS_PER_PIXEL*sizey_;
    int32_t centery = 0.5 + IMM_UNITS_PER_PIXEL*centery_;
    assert(sizex >= 0);
    assert(sizey >= 0);
    assert(centery >= 0);

    if ((sizex|centery) < (1 << 10) && (sizey) < (1 << 11))
    {
        int32_t btype;
        if (ibox_is_ptr(b))
        {
            btype = ibox_to_ptr(b)->get_type();
            delete ibox_to_ptr(b);
        }
        else
        {
            btype = iboximm_type(b);
        }
        b.imm =  uint64_t(1) + (uint64_t(sizey) << 1) + (uint64_t(centery) << 12) + (uint64_t(sizex) << 22)+ (uint64_t(btype) << 32);
    }
    else
    {
        if (ibox_is_ptr(b))
        {
            ibox_to_ptr(b)->sizex = sizex_;
            ibox_to_ptr(b)->sizey = sizey_;
            ibox_to_ptr(b)->centery = centery_;
        }
        else
        {
            b.ptr = new charbox(iboximm_type(b), sizex_, sizey_, centery_);
        }
    }

    return b;
}

ibox iboxchar_set_sizes(ibox b, double sizey_, double centery_)
{
    double sizex_ = ibox_is_ptr(b) ? ibox_to_ptr(b)->sizex : iboximm_sizex(b);
    int32_t sizex = 0.5 + IMM_UNITS_PER_PIXEL*sizex_;
    int32_t sizey = 0.5 + IMM_UNITS_PER_PIXEL*sizey_;
    int32_t centery = 0.5 + IMM_UNITS_PER_PIXEL*centery_;
    assert(sizex >= 0);
    assert(sizey >= 0);
    assert(centery >= 0);

    if ((sizex|centery) < (1 << 10) && (sizey) < (1 << 11))
    {
        int32_t btype;
        if (ibox_is_ptr(b))
        {
            btype = ibox_to_ptr(b)->get_type();
            delete ibox_to_ptr(b);
        }
        else
        {
            btype = iboximm_type(b);
        }
        b.imm =  uint64_t(1) + (uint64_t(sizey) << 1) + (uint64_t(centery) << 12) + (uint64_t(sizex) << 22)+ (uint64_t(btype) << 32);
    }
    else
    {
        if (ibox_is_ptr(b))
        {
            ibox_to_ptr(b)->sizex = sizex_;
            ibox_to_ptr(b)->sizey = sizey_;
            ibox_to_ptr(b)->centery = centery_;
        }
        else
        {
            b.ptr = new charbox(iboximm_type(b), sizex_, sizey_, centery_);
        }
    }

    return b;
}

rowbox * steal_rowbox(rowbox * row, int32_t a, int32_t b)
{
    size_t n = row->child.size();
    rowbox * newrow = new rowbox(n, a, b);
    for (int32_t i = 0; i < n; i++)
    {
        newrow->child[i].cibox = row->child[i].cibox;
        row->child[i].cibox = iboximm_make(0);
    }
    return newrow;
}

bool made_into_placeholder(rowbox * r)
{
    if (r->child.size() == 1)
    {
        r->child.insert(r->child.begin(), iboximm_make(CHAR_Placeholder));
        r->cursor_a = 0;
        r->cursor_b = 1;
        return true;
    }
    else
    {
        return false;
    }
}

int escape_seq_to_action(const char * s)
{
    int32_t c = esccode_to_char(s);
    if (c != 0)
    {
        return c;
    }
//  TODO: intt, sumt, ...
    return 0;
}
