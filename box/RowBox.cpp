#include <stack>
#include "boxes.h"
#include "ex_parse.h"
#include "notebook.h"
#include "box_convert.h"


rowbox::rowbox(const char * a, int32_t n)
     : boxbase(BNTYPE_ROW),
       style(gs.const_empty_list.copy()),
       options(gs.const_empty_list.copy())
{
    cursor_b = cursor_a = 0;
    int32_t i = 0;
    bool last_from_CR = false;
    while (i < n)
    {
        char16_t c;
        i += readonechar16(c, a + i);
        if (c == 10)
        {
            if (!last_from_CR)
            {
                child.push_back(iboxarrayelem(new nullbox()));
                cursor_b++;
            }
            last_from_CR = false;
        }
        else if (c == 13)
        {
            child.push_back(iboxarrayelem(new nullbox()));
            cursor_b++;
            last_from_CR = true;
        }
        else
        {
            child.push_back(iboxarrayelem(iboximm_make(c)));
            cursor_b++;
        }
    }
    child.push_back(iboxarrayelem(new nullbox()));
};


/* edit ******************************/

void rowbox::print(size_t indent, double offx, double offy)
{
    for (size_t i = 0; i < indent; i++)
        printf("    ");
    printf("(%f,%f) row: size(%f,%f:%f) cursor(%f,%f)\n",
             offx, offy, sizex, sizey, centery, cursor_a, cursor_b);

    for (auto& c: child)
    {
        if (ibox_is_imm(c.cibox))
        {
            for (size_t i = 0; i <= indent; i++)
                printf("    ");
            printf("(%f,%f) imm: (%d,%d,%d) size(%f,%f:%f)\n", c.offx, c.offy,
                 (iboximm_type(c.cibox)>>24)&255, (iboximm_type(c.cibox)>>16)&255, (iboximm_type(c.cibox)>>0)&65535,
                     iboximm_sizex(c.cibox), iboximm_sizey(c.cibox), iboximm_centery(c.cibox));
        }
        else
        {
            ibox_to_ptr(c.cibox)->print(indent + 1, c.offx, c.offy);
        }
    }
}

void rowbegin::print(size_t indent, double offx, double offy)
{
    for (size_t i = 0; i < indent; i++)
        printf("    ");

    printf("(%f,%f) rowbegin: size(%f,%f:%f)\n", offx, offy, sizex, sizey, sizey);
}

void rowend::print(size_t indent, double offx, double offy)
{
    for (size_t i = 0; i < indent; i++)
        printf("    ");

    printf("(%f,%f) rowend: size(%f,%f:%f)\n", offx, offy, sizex, sizey, sizey);
}


visitRet rowbegin::visit(visitArg m)
{
    return visitret_OK;
}

visitRet rowend::visit(visitArg m)
{
    return visitret_OK;
}

visitRet rowbox::visit(visitArg m)
{
    switch (m)
    {
        case visitarg_InvalidateAll:
        {
            flags &= ~(BNFLAG_MEASURED | BNFLAG_COLORED);
            for (auto i = child.begin(); i != child.end(); ++i)
            {
                if (ibox_is_ptr(i->cibox))
                    ibox_to_ptr(i->cibox)->visit(m);
            }
            return visitret_OK;
        }
        default:
        {
            assert(false);
            return visitret_OK;
        }
    }
}



boxbase * rowbegin::copy()
{
    rowbegin * r = new rowbegin();
    return r;
}

boxbase * rowend::copy()
{
    rowbegin * r = new rowbegin();
    return r;
}


boxbase * rowbox::copy()
{
    size_t n = child.size();
    rowbox * r = new rowbox(n, cursor_a, cursor_b);
    for (size_t i = 0; i < n; i++)
    {
        r->child[i].offx = child[i].offx;
        r->child[i].offy = child[i].offy;
        ibox b = child[i].cibox;
        if (ibox_is_ptr(b))
        {
            r->child[i].cibox.ptr = ibox_to_ptr(b)->copy();
        }
        else
        {
            r->child[i].cibox = b;
        }
    }
    return r;
}


void rowbegin::key_copy(boxbase*&b)
{
    return;
}

void rowend::key_copy(boxbase*&b)
{
    return;
}

void rowbox::key_copy(boxbase*&b)
{
    if (cursor_b >= child.size())
    {
        ibox_to_ptr(child[cursor_a].cibox)->key_copy(b);
    }
    else
    {
        int32_t x = std::min(cursor_a, cursor_b);
        int32_t y = std::max(cursor_a, cursor_b);
        if (x < y)
        {
            rowbox * r = new rowbox(y - x, 0, 0);
            for (int32_t i = x; i < y; i++)
                r->child[i - x].cibox = ibox_copy(child[i].cibox);
            assert(b == nullptr);
            b = r;
        }
    }
}

void rowbegin::key_paste(boxbase*&b)
{
    return;
}

void rowend::key_paste(boxbase*&b)
{
    return;
}

void rowbox::key_paste(boxbase*&b)
{
    flags &= ~(BNFLAG_MEASURED | BNFLAG_COLORED);

    if (cursor_b >= child.size())
    {
        ibox_to_ptr(child[cursor_a].cibox)->key_paste(b);
    }
    else if (b != nullptr && b->get_type() == BNTYPE_ROW)
    {
        delete_selection();
        rowbox* row = dynamic_cast<rowbox*>(b);
        size_t n = row->child.size();
        if (n > 0 && BNTYPE_NULLER == ibox_type(row->child[n - 1].cibox))
            n--;
        child.insert(child.begin() + cursor_a, n, iboxarrayelem(iboximm_make(0)));
        for (int32_t j = 0; j < n; j++)
            child[cursor_a + j].cibox = ibox_copy(row->child[j].cibox);
        cursor_b = cursor_a = cursor_a + n;
    }
}

void rowbox::delete_selection()
{
    assert(cursor_a < child.size());
    assert(cursor_b < child.size());

    int32_t left = std::min(cursor_a, cursor_b);
    int32_t right = std::max(cursor_a, cursor_b);

    for (auto i = child.begin() + left; i != child.begin() + right; ++i)
        if (ibox_is_ptr((*i).cibox))
            delete ibox_to_ptr((*i).cibox);
    child.erase(child.begin() + left, child.begin() + right);
    cursor_b = cursor_a = left;
}

void rowbegin::insert_char(int32_t c)
{
    assert(false);
}

void rowend::insert_char(int32_t c)
{
    assert(false);
}

void rowbox::insert_char(int32_t c)
{
    flags &= ~(BNFLAG_MEASURED | BNFLAG_COLORED);

    if (cursor_b >= child.size())
    {
        assert(cursor_a < child.size());
        ibox_to_ptr(child[cursor_a].cibox)->insert_char(c);
        return;
    }

    if (c < 0)
    {
        if ((c == -11 || c == -12) && cursor_a != cursor_b)
        {
            int32_t left = std::min(cursor_a, cursor_b);
            int32_t right = std::max(cursor_a, cursor_b);
            rowbox * newrow = new rowbox(1 + right - left, 0, 0);
            for (int32_t i = left; i < right; ++i)
            {
                newrow->child[i-left].cibox = child[i].cibox;
                child[i].cibox = iboximm_make(0);

            }
            newrow->child[right-left].cibox.ptr = new nullbox();
            uint32_t angle = 1<<24;
            if (c == -11)
                angle = -angle;
            child[left].cibox.ptr = new rotationbox(newrow, angle, angle);
            child.erase(child.begin() + left + 1, child.begin() + right);
            cursor_a = left;
            cursor_b = child.size();
        }
        return;
    }

    delete_selection();
    assert(cursor_a < child.size());
    child.insert(child.begin() + cursor_a, iboxarrayelem());
    child[cursor_a].cibox = iboximm_make(c);
    cursor_a++;
    cursor_b = cursor_a;

    int32_t r, l, cc;
    char t[40];

    // look for \[ to the left
    for (l = cursor_a - 1; l >= 0; l--)
    {
        if (!ibox_is_char(child[l].cibox))
            goto alias_scan;
        if (ibox_is_char(child[l].cibox, '\\'))
            break;
    }
    if (l < 0)
        goto alias_scan;
    if (!(l+2 < child.size() && ibox_is_char(child[l+1].cibox, '[')))
        goto alias_scan;
    if (l > 0 && ibox_is_char(child[l-1].cibox, '['))
        goto alias_scan;

    // look for ] to the right
    for (r = cursor_a - 1; r+1 < child.size(); r++)
    {
        if (!ibox_is_char(child[r].cibox))
            goto alias_scan;
        if (ibox_is_char(child[r].cibox, ']'))
            break;
    }
    if (r+1 >= child.size())
        goto alias_scan;

    if (r-l > 30)
        goto alias_scan;

    for (int32_t i=l; i<=r; i++)
        t[i-l] = ibox_type(child[i].cibox)&65535;
    t[r+1-l] = 0;
    cc = escapedname_to_char(t) && 65535;
    if (cc == 0)
        goto alias_scan;

    for (auto i = child.begin() + l; i != child.begin() + r+1; ++i)
        if (ibox_is_ptr((*i).cibox))
            delete ibox_to_ptr((*i).cibox);
    child.erase(child.begin() + l+1, child.begin() + r+1);
    child[l].cibox = iboximm_make(cc);
    cursor_b = cursor_a = l + 1;
    return;

alias_scan:

    // look for \[AliasDelimiter] to the right
    for (r = cursor_a - 1; r+1 < child.size(); r++)
    {
        if (!ibox_is_char(child[r].cibox))
            return;
        if (ibox_is_char(child[r].cibox, CHAR_AliasDelimiter))
            break;
    }
    if (r+1 >= child.size())
        return;

    // look for \[AliasDelimiter] to the left
    for (l = std::min(cursor_a - 1, r-1); l >= 0; l--)
    {
        if (!ibox_is_char(child[l].cibox))
            return;
        if (ibox_is_char(child[l].cibox, CHAR_AliasDelimiter))
            break;
    }
    if (l<0)
        return;

    if (r-l > 30)
        return;

    for (int32_t i=l+1; i<r; i++)
        t[i-(l+1)] = ibox_type(child[i].cibox)&65535;
    t[r-(l+1)] = 0;
    cc = escape_seq_to_action(t);
    if (cc <= 0)
        return;

    for (auto i = child.begin() + l; i != child.begin() + r+1; ++i)
        if (ibox_is_ptr((*i).cibox))
            delete ibox_to_ptr((*i).cibox);
    child.erase(child.begin() + l+1, child.begin() + r+1);
    child[l].cibox = iboximm_make(cc);
    cursor_b = cursor_a = l + 1;
    return;
}


moveRet rowbegin::move(boxbase*&b, moveArg m)
{
    assert(false);
    return moveret_OK;
}

moveRet rowend::move(boxbase*&b, moveArg m)
{
    assert(false);
    return moveret_OK;
}

moveRet rowbox::move(boxbase*&b, moveArg m)
{
    moveRet r;
    switch (m)
    {
        case movearg_Left:
        {
            if (cursor_b >= child.size())
            {
                assert(cursor_a < child.size());
                r = ibox_to_ptr(child[cursor_a].cibox)->move(b, movearg_Left);
                assert(r == moveret_OK || r == moveret_End);
                if (r == moveret_End)
                    cursor_b = cursor_a;
                return moveret_OK;
            }
            else
            {
                assert(cursor_a < child.size());
                if (cursor_a > 0)
                {
                    cursor_a--;
                    r = moveret_End;
                    if (ibox_is_ptr(child[cursor_a].cibox))
                        r = ibox_to_ptr(child[cursor_a].cibox)->move(b, movearg_Last);
                    assert(r == moveret_OK || r == moveret_End);
                    if (r == moveret_End)
                        cursor_b = cursor_a;
                    else
                        cursor_b = child.size();
                    return moveret_OK;
                }
                else
                {
                    cursor_a = cursor_b = 0;
                    return moveret_End;
                }
            }
        }
        case movearg_ShiftLeft:
        {
            if (cursor_b >= child.size())
            {
                assert(cursor_a < child.size());
                r = ibox_to_ptr(child[cursor_a].cibox)->move(b, movearg_ShiftLeft);
                assert(r == moveret_OK || r == moveret_End);
                if (r == moveret_End)
                    cursor_b = std::min(cursor_a + 1, int32_t(child.size() - 1));
                return moveret_OK;
            }
            else
            {
                assert(cursor_a < child.size());
                if (cursor_a > 0)
                {
                    cursor_a--;
                    return moveret_OK;
                }
                else
                {
                    return moveret_End;
                }
            }
        }
        case movearg_Right:
        {
            if (cursor_b >= child.size())
            {
                assert(cursor_a < child.size());
                r = ibox_to_ptr(child[cursor_a].cibox)->move(b, movearg_Right);
            }
            else
            {
                assert(cursor_a < child.size());
                r = moveret_End;
                if (ibox_is_ptr(child[cursor_a].cibox))
                    r = ibox_to_ptr(child[cursor_a].cibox)->move(b, movearg_First);
            }
            assert(r == moveret_OK || r == moveret_End);
            if (r == moveret_End)
            {
                if (cursor_a + 1 < child.size())
                {
                    cursor_a++;
                    cursor_b = cursor_a;
                    return moveret_OK;
                }
                else
                {
                    cursor_b = cursor_a;
                    return moveret_End;
                }
            }
            else
            {
                cursor_b = child.size();
                return moveret_OK;
            }
        }
        case movearg_ShiftRight:
        {
            if (cursor_b >= child.size())
            {
                assert(cursor_a < child.size());
                r = ibox_to_ptr(child[cursor_a].cibox)->move(b, movearg_ShiftRight);
                assert(r == moveret_OK || r == moveret_End);
                if (r == moveret_OK)
                    return r;
                cursor_b = cursor_a;
            }
            if (cursor_a + 1 < child.size())
            {
                cursor_a++;
                return moveret_OK;
            }
            else
            {
                return moveret_End;
            }
        }
        case movearg_ShiftUp:
        case movearg_Up:
        {
            if (cursor_b >= child.size())
            {
                assert(cursor_a < child.size());
                r = ibox_to_ptr(child[cursor_a].cibox)->move(b, m);
                assert(r == moveret_OK || r == moveret_End);
                if (r == moveret_OK)
                    return r;
                cursor_b = cursor_a;
            }
            int32_t old_b = cursor_b;
            int32_t cursorx = child[cursor_a].offx;
            int32_t leftx = cursorx;
            int32_t rightx = cursorx;
            int line_count = 0;
            for (int32_t t = cursor_a - 1; t >= 0; t--)
            {
                rightx = leftx;
                leftx = child[t].offx;
                if (leftx >= rightx)
                {
                    line_count++;
                    if (leftx <= cursorx)
                    {
                        cursor_b = cursor_a = t;
                        if (m == movearg_ShiftUp)
                            cursor_b = old_b;
                        return moveret_OK;
                    }
                }
                else if (line_count < 1)
                {
                    continue;
                }
                if (leftx <= cursorx && cursorx <= rightx)
                {
                    if (cursorx - leftx > rightx - cursorx)
                        t++;
                    cursor_b = cursor_a = t;
                    if (m == movearg_ShiftUp)
                        cursor_b = old_b;
                    return moveret_OK;
                }
            }

            cursor_b = cursor_a;
            if (m == movearg_ShiftUp)
                cursor_b = old_b;
            return moveret_End;
        }
        case movearg_ShiftDown:
        case movearg_Down:
        {
            if (cursor_b >= child.size())
            {
                assert(cursor_a < child.size());
                r = ibox_to_ptr(child[cursor_a].cibox)->move(b, m);
                assert(r == moveret_OK || r == moveret_End);
                if (r == moveret_OK)
                    return r;
                cursor_b = cursor_a;
            }
            int32_t old_b = cursor_b;
            int32_t cursorx = child[cursor_a].offx;
            int32_t leftx = cursorx;
            int32_t rightx = cursorx;
            int line_count = 0;
            for (int32_t t = cursor_a + 1; t < child.size(); t++)
            {
                leftx = rightx;
                rightx = child[t].offx;
                if (leftx >= rightx)
                {
                    line_count++;
                    if (line_count > 1)
                    {
                        cursor_b = cursor_a = t - 1;
                        if (m == movearg_ShiftDown)
                            cursor_b = old_b;
                        return moveret_OK;
                    }
                }
                if (line_count < 1)
                    continue;
                if (leftx <= cursorx && cursorx <= rightx)
                {
                    if (cursorx - leftx <= rightx - cursorx && t - 1 > cursor_a)
                        t--;
                    cursor_b = cursor_a = t;
                    if (m == movearg_ShiftDown)
                        cursor_b = old_b;
                    return moveret_OK;
                }
            }

            if (line_count > 0)
            {
                cursor_b = cursor_a = child.size() - 1;
                if (m == movearg_ShiftDown)
                    cursor_b = old_b;
                return moveret_OK;
            }
            else
            {
                cursor_b = cursor_a;
                if (m == movearg_ShiftDown)
                    cursor_b = old_b;
                return moveret_End;
            }
        }
        case movearg_Last:
        {
            assert(child.size() > 0);
            cursor_b = cursor_a = child.size() - 1;
            if (cursor_a > 0 &&
                ibox_is_ptr(child[cursor_a - 1].cibox) &&
                ibox_to_ptr(child[cursor_a - 1].cibox)->get_type() == BNTYPE_MONO)
            {
                cursor_b = child.size();
                cursor_a--;
                ibox_to_ptr(child[cursor_a].cibox)->move(b, movearg_Last);
            }
            else
            {
                select_placeholder_if_possible();
            }
            return moveret_OK;
        }
        case movearg_First:
        {
            assert(child.size() > 0);
            cursor_b = cursor_a = 0;
            if (cursor_a < child.size() - 1 &&
                ibox_is_ptr(child[cursor_a].cibox) &&
                ibox_to_ptr(child[cursor_a].cibox)->get_type() == BNTYPE_MONO)
            {
                cursor_b = child.size();
                ibox_to_ptr(child[cursor_a].cibox)->move(b, movearg_First);
            }
            else
            {
                select_placeholder_if_possible();
            }
            return moveret_OK;
        }
        case movearg_Home:
        {
            if (cursor_b >= child.size())
            {
                assert(cursor_a < child.size());
                r = ibox_to_ptr(child[cursor_a].cibox)->move(b, m);
                if (r == moveret_End)
                    cursor_b = cursor_a;
                return r;
            }
            else
            {
                r = moveret_End;
                while (cursor_a > 0 && ibox_type(child[cursor_a - 1].cibox) != BNTYPE_NULLER)
                {
                    cursor_a--;
                    cursor_b = cursor_a;
                    r = moveret_OK;
                }
                return r;
            }
        }
        case movearg_End:
        {
            if (cursor_b >= child.size())
            {
                assert(cursor_a < child.size());
                r = ibox_to_ptr(child[cursor_a].cibox)->move(b, m);
                if (r == moveret_End)
                    cursor_b = cursor_a = cursor_a + 1;
                return r;
            }
            else
            {
                r = moveret_End;
                while (cursor_a + 1 < child.size() && ibox_type(child[cursor_a].cibox) != BNTYPE_NULLER)
                {
                    cursor_b = cursor_a = cursor_a + 1;
                    r = moveret_OK;
                }
                return r;
            }
        }
        case movearg_CtrlSpace:
        {
            if (cursor_b >= child.size())
            {
                assert(cursor_a < child.size());
                r = ibox_to_ptr(child[cursor_a].cibox)->move(b, m);
                if (r == moveret_End)
                    cursor_b = cursor_a = cursor_a + 1;
                return moveret_OK;
            }
            else
            {
                return moveret_End;
            }
        }
        case movearg_Tab:
        {
            flags &= ~(BNFLAG_MEASURED | BNFLAG_COLORED);
            if (cursor_b >= child.size())
            {
                assert(cursor_a < child.size());
                r = ibox_to_ptr(child[cursor_a].cibox)->move(b, m);
                if (r == moveret_End)
                    cursor_b = cursor_a = cursor_a + 1;
                return r;
            }
            else
            {
                delete_selection();
                for (int i = 0; i < 4; i++)
                {
                    assert(cursor_a < child.size());
                    child.insert(child.begin() + cursor_a, iboxarrayelem());
                    child[cursor_a].cibox = iboximm_make(' ');
                    cursor_b = cursor_a = cursor_a + 1;
                }
                return moveret_OK;
            }
        }
        case movearg_Switch:
        {
            if (cursor_b >= child.size())
            {
                flags &= ~(BNFLAG_MEASURED | BNFLAG_COLORED);
                assert(cursor_a < child.size());
                r = ibox_to_ptr(child[cursor_a].cibox)->move(b, movearg_Switch);
                if (r == moveret_OK)
                {
                    return r;
                }
                else if (r == moveret_Replace)
                {
                    replacechild(cursor_a, b);
                    b = nullptr;
                    return moveret_OK;
                }
                else
                {
                    assert(false);
                    return r;
                }
            }
            else
            {
                return moveret_OK;
            }
        }
        default:
        {
            assert(false);
            return moveret_OK;
        }
    }
}

void rowbox::select_placeholder_if_possible()
{
    if (child.size() == 2 && ibox_is_char(child[0].cibox, CHAR_Placeholder))
    {
        cursor_b = cursor_a == 0 ? 1 : 0;
    }
}

void rowbox::select_prev_if_possible()
{
    int32_t npi, i;
    int32_t pi = cursor_a;

    npi = pi;
    while (npi > 0)
    {
        int32_t mtype, ptype;
        if (!ibox_is_char(child[npi - 1].cibox))
        {
            if (pi == npi && ibox_is_node(child[npi - 1].cibox))
                npi--;
            break;
        }
        ptype = ibox_type(child[npi - 1].cibox)&65535;
        if (ptype == ']' || ptype == '}' || ptype == ')' || ptype == CHAR_RightDoubleBracket)
        {
            std::stack<int32_t> estack;
            estack.push(ptype);
            while (--npi > 0 && estack.size() > 0)
            {
                if (!ibox_is_char(child[npi - 1].cibox)) {
                    continue;
                }
                ptype = ibox_type(child[npi - 1].cibox)&65535;
                if (ptype == ']' || ptype == '}' || ptype == ')' || ptype == CHAR_RightDoubleBracket) {
                    mtype = 0;
                    estack.push(ptype);
                } else if (ptype == '[') {
                    mtype = ']';
                } else if (ptype == CHAR_LeftDoubleBracket) {
                    mtype = CHAR_RightDoubleBracket;
                } else if (ptype == '(') {
                    mtype = ')';
                } else if (ptype == '{') {
                    mtype = '}';
                } else {
                    continue;
                }
                if (estack.top() == mtype) {
                    estack.pop();
                } else {
                    estack.push(ptype);
                }                
            }
        }
        else if (isletterchar(ptype) || ('0' <= ptype && ptype <= '9'))
        {
            --npi;
        }
        else if (pi == npi && (ptype == CHAR_Sum || ptype == CHAR_Product))
        {
            --npi;
            break;
        }
        else
        {
            break;
        }
    }
    if (npi < pi)
    {
        cursor_b = npi;
    }
}


insertRet rowbegin::insert(boxbase*&b, insertArg m)
{
    assert(false);
    return insertret_Done;
}

insertRet rowend::insert(boxbase*&b, insertArg m)
{
    assert(false);
    return insertret_Done;
}

insertRet rowbox::insert(boxbase*&b, insertArg m)
{
//printf("rowbox::insert called m = %d\n", m);

    flags &= ~(BNFLAG_MEASURED | BNFLAG_COLORED);

    insertRet r;
    switch (m)
    {
        case insertarg_Fraction:
        {
            if (cursor_b >= child.size())
            {
                return ibox_to_ptr(child[cursor_a].cibox)->insert(b, m);
            }
            else if (cursor_a == cursor_b)
            {
                select_prev_if_possible();
            }
            if (cursor_a == cursor_b)
            {
                rowbox* newrow1 = new rowbox(2, 0,1);
                newrow1->child[0].cibox = iboximm_make(CHAR_Placeholder);
                newrow1->child[1].cibox.ptr = new nullbox();
                rowbox* newrow2 = new rowbox(2, 0,0);
                newrow2->child[0].cibox = iboximm_make(CHAR_Placeholder);
                newrow2->child[1].cibox.ptr = new nullbox();
                fractionbox* newfraction = new fractionbox(newrow1, newrow2, 0);
                child.insert(child.begin() + cursor_a, iboxarrayelem(newfraction));
                cursor_b = child.size();
                return insertret_Done;
            }
            else
            {
                int32_t left = std::min(cursor_a, cursor_b);
                int32_t right = std::max(cursor_a, cursor_b);
                rowbox* newrow1 = new rowbox(right - left + 1, 0,0);
                for (int32_t i = left; i < right; i++)
                    newrow1->child[i - left].cibox = child[i].cibox;
                newrow1->child[right - left].cibox.ptr = new nullbox();
                child.erase(child.begin() + left, child.begin() + right);
                rowbox* newrow2 = new rowbox(2, 0,1);
                newrow2->child[0].cibox = iboximm_make(CHAR_Placeholder);
                newrow2->child[1].cibox.ptr = new nullbox();
                fractionbox* newfraction = new fractionbox(newrow1, newrow2, 1);                
                child.insert(child.begin() + left, iboxarrayelem(newfraction));
                cursor_a = left;
                cursor_b = child.size();
                return insertret_Done;
            }
        }
        case insertarg_Sqrt:
        {
            if (cursor_b >= child.size())
            {
                return ibox_to_ptr(child[cursor_a].cibox)->insert(b, m);
            }
            else if (cursor_a == cursor_b)
            {
                rowbox* newrow1 = new rowbox(2, 0,1);
                newrow1->child[0].cibox = iboximm_make(CHAR_Placeholder);
                newrow1->child[1].cibox.ptr = new nullbox();
                sqrtbox* newsqrt = new sqrtbox(newrow1);
                child.insert(child.begin() + cursor_a, iboxarrayelem(newsqrt));
                cursor_b = child.size();
                return insertret_Done;
            }
            else
            {
                int32_t left = std::min(cursor_a, cursor_b);
                int32_t right = std::max(cursor_a, cursor_b);
                rowbox* newrow1 = new rowbox(right - left + 1, 0, (right - left == 1 && ibox_is_char(child[0].cibox)) ? 1 : 0);
                for (int32_t i = left; i < right; i++)
                    newrow1->child[i - left].cibox = child[i].cibox;
                newrow1->child[right - left].cibox.ptr = new nullbox();
                child.erase(child.begin() + left, child.begin() + right);
                sqrtbox* newsqrt = new sqrtbox(newrow1);                
                child.insert(child.begin() + left, iboxarrayelem(newsqrt));
                cursor_a = left;
                cursor_b = child.size();
                return insertret_Done;
            }
        }
        case insertarg_Newline:
        {
            if (cursor_b >= child.size())
            {
                return ibox_to_ptr(child[cursor_a].cibox)->insert(b, m);
            }
            delete_selection();
            child.insert(child.begin() + cursor_a, iboxarrayelem(new nullbox()));
            cursor_b = cursor_a = cursor_a + 1;
            return insertret_Done;
        }
        case insertarg_Subscript:
        {
            if (cursor_b >= child.size())
            {
                return ibox_to_ptr(child[cursor_a].cibox)->insert(b, m);
            }
            else if (cursor_a == cursor_b)
            {
                rowbox* newrow1 = new rowbox(2, 0,1);
                newrow1->child[0].cibox = iboximm_make(CHAR_Placeholder);
                newrow1->child[1].cibox.ptr = new nullbox();
                subscriptbox* newsubscript = new subscriptbox(newrow1);
                child.insert(child.begin() + cursor_a, iboxarrayelem(newsubscript));
                cursor_b = child.size();
                return insertret_Done;
            }
            else
            {
                int32_t left = std::min(cursor_a, cursor_b);
                int32_t right = std::max(cursor_a, cursor_b);
                rowbox* newrow1 = new rowbox(right - left + 1, 0, (right - left == 1 && ibox_is_char(child[0].cibox)) ? 1 : 0);
                for (int32_t i = left; i < right; i++)
                    newrow1->child[i - left].cibox = child[i].cibox;
                newrow1->child[right - left].cibox.ptr = new nullbox();
                child.erase(child.begin() + left, child.begin() + right);
                subscriptbox* newsubscript = new subscriptbox(newrow1);
                child.insert(child.begin() + left, iboxarrayelem(newsubscript));
                cursor_a = left;
                cursor_b = child.size();
                return insertret_Done;
            }
        }
        case insertarg_Superscript:
        {
            if (cursor_b >= child.size())
            {
                return ibox_to_ptr(child[cursor_a].cibox)->insert(b, m);
            }
            else if (cursor_a == cursor_b)
            {
                rowbox* newrow1 = new rowbox(2, 0,1);
                newrow1->child[0].cibox = iboximm_make(CHAR_Placeholder);
                newrow1->child[1].cibox.ptr = new nullbox();
                superscriptbox* newsuperscript = new superscriptbox(newrow1);
                child.insert(child.begin() + cursor_a, iboxarrayelem(newsuperscript));
                cursor_b = child.size();
                return insertret_Done;
            }
            else
            {
                int32_t left = std::min(cursor_a, cursor_b);
                int32_t right = std::max(cursor_a, cursor_b);
                rowbox* newrow1 = new rowbox(right - left + 1, 0, (right - left == 1 && ibox_is_char(child[0].cibox)) ? 1 : 0);
                for (int32_t i = left; i < right; i++)
                    newrow1->child[i - left].cibox = child[i].cibox;
                newrow1->child[right - left].cibox.ptr = new nullbox();
                child.erase(child.begin() + left, child.begin() + right);
                superscriptbox* newsuperscript = new superscriptbox(newrow1);
                child.insert(child.begin() + left, iboxarrayelem(newsuperscript));
                cursor_a = left;
                cursor_b = child.size();
                return insertret_Done;
            }
        }
        case insertarg_Underscript:
        {
            if (cursor_b >= child.size())
            {
                return ibox_to_ptr(child[cursor_a].cibox)->insert(b, m);
            }
            else if (cursor_a == cursor_b)
            {
                select_prev_if_possible();
            }
            if (cursor_a == cursor_b)
            {
                rowbox* newrow1 = new rowbox(2, 0,1);
                newrow1->child[0].cibox = iboximm_make(CHAR_Placeholder);
                newrow1->child[1].cibox.ptr = new nullbox();
                rowbox* newrow2 = new rowbox(2, 0,0);
                newrow2->child[0].cibox = iboximm_make(CHAR_Placeholder);
                newrow2->child[1].cibox.ptr = new nullbox();
                underscriptbox* newunderscript = new underscriptbox(newrow1, newrow2, 0);
                child.insert(child.begin() + cursor_a, iboxarrayelem(newunderscript));
                cursor_b = child.size();
                return insertret_Done;
            }
            else
            {
                int32_t left = std::min(cursor_a, cursor_b);
                int32_t right = std::max(cursor_a, cursor_b);
                rowbox* newrow1 = new rowbox(right - left + 1, 0,0);
                for (int32_t i = left; i < right; i++)
                    newrow1->child[i - left].cibox = child[i].cibox;
                newrow1->child[right - left].cibox.ptr = new nullbox();
                child.erase(child.begin() + left, child.begin() + right);
                rowbox* newrow2 = new rowbox(2, 0,1);
                newrow2->child[0].cibox = iboximm_make(CHAR_Placeholder);
                newrow2->child[1].cibox.ptr = new nullbox();
                underscriptbox* newunderscript = new underscriptbox(newrow1, newrow2, 1);                
                child.insert(child.begin() + left, iboxarrayelem(newunderscript));
                cursor_a = left;
                cursor_b = child.size();
                return insertret_Done;
            }
        }
        case insertarg_Overscript:
        {
            if (cursor_b >= child.size())
            {
                return ibox_to_ptr(child[cursor_a].cibox)->insert(b, m);
            }
            else if (cursor_a == cursor_b)
            {
                select_prev_if_possible();
            }
            if (cursor_a == cursor_b)
            {
                rowbox* newrow1 = new rowbox(2, 0,1);
                newrow1->child[0].cibox = iboximm_make(CHAR_Placeholder);
                newrow1->child[1].cibox.ptr = new nullbox();
                rowbox* newrow2 = new rowbox(2, 0,0);
                newrow2->child[0].cibox = iboximm_make(CHAR_Placeholder);
                newrow2->child[1].cibox.ptr = new nullbox();
                overscriptbox* newoverscript = new overscriptbox(newrow1, newrow2, 0);
                child.insert(child.begin() + cursor_a, iboxarrayelem(newoverscript));
                cursor_b = child.size();
                return insertret_Done;
            }
            else
            {
                int32_t left = std::min(cursor_a, cursor_b);
                int32_t right = std::max(cursor_a, cursor_b);
                rowbox* newrow1 = new rowbox(right - left + 1, 0,0);
                for (int32_t i = left; i < right; i++)
                    newrow1->child[i - left].cibox = child[i].cibox;
                newrow1->child[right - left].cibox.ptr = new nullbox();
                child.erase(child.begin() + left, child.begin() + right);
                rowbox* newrow2 = new rowbox(2, 0,1);
                newrow2->child[0].cibox = iboximm_make(CHAR_Placeholder);
                newrow2->child[1].cibox.ptr = new nullbox();
                overscriptbox* newoverscript = new overscriptbox(newrow1, newrow2, 1);                
                child.insert(child.begin() + left, iboxarrayelem(newoverscript));
                cursor_a = left;
                cursor_b = child.size();
                return insertret_Done;
            }
        }
        case insertarg_GridRow:
        {
            if (cursor_b >= child.size())
            {
                return ibox_to_ptr(child[cursor_a].cibox)->insert(b, m);
            }
            else if (cursor_a == cursor_b)
            {
                select_prev_if_possible();
            }
            if (cursor_a == cursor_b)
            {
                rowbox* newrow1 = new rowbox(2, 0,1);
                newrow1->child[0].cibox = iboximm_make(CHAR_Placeholder);
                newrow1->child[1].cibox.ptr = new nullbox();
                rowbox* newrow2 = new rowbox(2, 0,0);
                newrow2->child[0].cibox = iboximm_make(CHAR_Placeholder);
                newrow2->child[1].cibox.ptr = new nullbox();
                gridbox* newgrid = new gridbox(2,1, 1,0);
                newgrid->array[0][0].cbox = newrow1;
                newgrid->array[1][0].cbox = newrow2;
                child.insert(child.begin() + cursor_a, iboxarrayelem(newgrid));
                cursor_b = child.size();
                return insertret_Done;
            }
            else
            {
                int32_t left = std::min(cursor_a, cursor_b);
                int32_t right = std::max(cursor_a, cursor_b);
                rowbox* newrow1 = new rowbox(right - left + 1, 0,0);
                for (int32_t i = left; i < right; i++)
                    newrow1->child[i - left].cibox = child[i].cibox;
                newrow1->child[right - left].cibox.ptr = new nullbox();
                child.erase(child.begin() + left, child.begin() + right);
                rowbox* newrow2 = new rowbox(2, 0,1);
                newrow2->child[0].cibox = iboximm_make(CHAR_Placeholder);
                newrow2->child[1].cibox.ptr = new nullbox();
                gridbox* newgrid = new gridbox(2,1, 1,0);
                newgrid->array[0][0].cbox = newrow1;
                newgrid->array[1][0].cbox = newrow2;
                child.insert(child.begin() + left, iboxarrayelem(newgrid));
                cursor_a = left;
                cursor_b = child.size();
                return insertret_Done;
            }
        }
        case insertarg_GridCol:
        {
            if (cursor_b >= child.size())
            {
                return ibox_to_ptr(child[cursor_a].cibox)->insert(b, m);
            }
            else if (cursor_a == cursor_b)
            {
                select_prev_if_possible();
            }
            if (cursor_a == cursor_b)
            {
                rowbox* newrow1 = new rowbox(2, 0,1);
                newrow1->child[0].cibox = iboximm_make(CHAR_Placeholder);
                newrow1->child[1].cibox.ptr = new nullbox();
                rowbox* newrow2 = new rowbox(2, 0,0);
                newrow2->child[0].cibox = iboximm_make(CHAR_Placeholder);
                newrow2->child[1].cibox.ptr = new nullbox();
                gridbox* newgrid = new gridbox(1,2, 0,1);
                newgrid->array[0][0].cbox = newrow1;
                newgrid->array[0][1].cbox = newrow2;
                child.insert(child.begin() + cursor_a, iboxarrayelem(newgrid));
                cursor_b = child.size();
                return insertret_Done;
            }
            else
            {
                int32_t left = std::min(cursor_a, cursor_b);
                int32_t right = std::max(cursor_a, cursor_b);
                rowbox* newrow1 = new rowbox(right - left + 1, 0,0);
                for (int32_t i = left; i < right; i++)
                    newrow1->child[i - left].cibox = child[i].cibox;
                newrow1->child[right - left].cibox.ptr = new nullbox();
                child.erase(child.begin() + left, child.begin() + right);
                rowbox* newrow2 = new rowbox(2, 0,1);
                newrow2->child[0].cibox = iboximm_make(CHAR_Placeholder);
                newrow2->child[1].cibox.ptr = new nullbox();
                gridbox* newgrid = new gridbox(1,2, 0,1);
                newgrid->array[0][0].cbox = newrow1;
                newgrid->array[0][1].cbox = newrow2;
                child.insert(child.begin() + left, iboxarrayelem(newgrid));
                cursor_a = left;
                cursor_b = child.size();
                return insertret_Done;
            }
        }
        case insertarg_Text:
        {
            if (cursor_b >= child.size())
            {
                return ibox_to_ptr(child[cursor_a].cibox)->insert(b, m);
            }
            delete_selection();
            child.insert(child.begin() + cursor_a, iboxarrayelem(new monobox(1,0,0, 0,0)));
            cursor_b = child.size();
            return insertret_Done;
        }
        case insertarg_Graphics3D:
        {
            if (cursor_b >= child.size())
            {
                return ibox_to_ptr(child[cursor_a].cibox)->insert(b, m);
            }
            delete_selection();
            child.insert(child.begin() + cursor_a, iboxarrayelem(new graphics3dbox()));
            cursor_b = cursor_a = cursor_a + 1;
            return insertret_Done;
        }
        default:
        {
            assert(false);
            return insertret_Done;
        }
    }
}

removeRet rowbegin::remove(boxbase*&b, removeArg m)
{
    assert(false);
    return removeret_OK;
}

removeRet rowend::remove(boxbase*&b, removeArg m)
{
    assert(false);
    return removeret_OK;
}

removeRet rowbox::remove(boxbase*&b, removeArg m)
{
//std::cout << "rowbox::remove called" << std::endl;
    assert(b == nullptr);
    removeRet r;
    moveRet s;

    flags &= ~(BNFLAG_MEASURED | BNFLAG_COLORED);

    if (cursor_b >= child.size())
    {
        assert(cursor_a < child.size());
        r = ibox_to_ptr(child[cursor_a].cibox)->remove(b, m);
        assert(r == removeret_OK || r == removeret_Replace || r == removeret_End);
        if (r == removeret_Replace)
        {
            delete ibox_to_ptr(child[cursor_a].cibox);
            if (b == nullptr)
            {
                child.erase(child.begin() + cursor_a);
                cursor_b = cursor_a;
                return removeret_OK;
            }
            else if (b->get_type() == BNTYPE_ROW)
            {
                child.erase(child.begin() + cursor_a);
                cursor_b = cursor_a;
                rowbox* row = dynamic_cast<rowbox*>(b);
                size_t n = row->child.size();
                if (n > 0 && BNTYPE_NULLER == ibox_type(row->child[n - 1].cibox))
                    n--;
                child.insert(child.begin() + cursor_a, n, iboxarrayelem(iboximm_make(0)));
                for (int32_t j = 0; j < n; j++)
                {
                    child[cursor_a + j].cibox = row->child[j].cibox;
                    row->child[j].cibox = iboximm_make(0);
                }
                if (m == removearg_Left)
                    cursor_b = cursor_a = cursor_a + n;
                delete b;
                b = nullptr;
                return removeret_OK;
            }
            else
            {
                child[cursor_a].cibox.ptr = b;
                b = nullptr;
                return removeret_OK;
            }
        }
        else if (r == removeret_End)
        {
            if (m == removearg_Right)
                cursor_a++;
            cursor_b = cursor_a;
            return removeret_OK;
        }
        else
        {
            return removeret_OK;
        }
    }
    else if (cursor_a == cursor_b)
    {
        assert(cursor_a < child.size());
        switch (m)
        {
            case removearg_Left:
            {
                if (cursor_a <= 0)
                {
                    return removeret_End;
                }
                else if (ibox_is_char(child[cursor_a - 1].cibox))
                {
                    if (ibox_is_ptr(child[cursor_a - 1].cibox))
                        delete ibox_to_ptr(child[cursor_a - 1].cibox);
                    child.erase(child.begin() + cursor_a - 1);
                    cursor_b = cursor_a = cursor_a - 1;
                    return removeret_OK;
                }
                else
                {
                    cursor_a--;
                    s = ibox_to_ptr(child[cursor_a].cibox)->move(b, movearg_Last);
                    assert(s == moveret_OK || s == moveret_End);
                    if (s == moveret_End)
                    {
                        delete ibox_to_ptr(child[cursor_a].cibox);
                        child.erase(child.begin() + cursor_a);
                        cursor_b = cursor_a;
                    }
                    else
                    {
                        cursor_b = child.size();
                    }
                    return removeret_OK;
                }
            }
            case removearg_Right:
            {
                if (cursor_a + 1 >= child.size())
                {
                    return removeret_End;
                }
                else if (ibox_is_char(child[cursor_a].cibox))
                {
                    if (ibox_is_ptr(child[cursor_a].cibox))
                        delete ibox_to_ptr(child[cursor_a].cibox);
                    child.erase(child.begin() + cursor_a);
                    return removeret_OK;
                }
                else
                {
                    s = ibox_to_ptr(child[cursor_a].cibox)->move(b, movearg_First);
                    assert(s == moveret_OK || s == moveret_End);
                    if (s == moveret_End)
                    {
                        delete ibox_to_ptr(child[cursor_a].cibox);
                        child.erase(child.begin() + cursor_a);
                        return removeret_OK;
                    }
                    else
                    {
                        cursor_b = child.size();
                        return removeret_OK;
                    }
                }
            }
            default:
            {
                assert(false);
                return removeret_Bad;
            }
        }
    }
    else
    {
        delete_selection();
        return removeret_OK;
    }
}

ex rowbegin::get_ex()
{
    assert(false);
    return nullptr;
}

ex rowend::get_ex()
{
    assert(false);
    return nullptr;
}

static ex finish_row(std::vector<wex>&v)
{
    if (v.empty())
    {
        return emake_str("", 0);
    }
    else if (v.size() == 1 && eis_str(v[0].get()))
    {
        return ecopy(v[0].get());
    }
    else
    {
        ex t = emake_node(gs.sym_sList.copy(), v);
        return emake_node(gs.sym_sRowBox.copy(), t);
    }
}

ex rowbox::get_ex()
{
    std::vector<std::vector<wex>> w;
    w.push_back(std::vector<wex>());

    int32_t n = child.size();
    if (n > 0 && !ibox_is_char(child[n-1].cibox) &&
        ibox_to_ptr(child[n-1].cibox)->get_type() == BNTYPE_NULLER)
    {
        n--;
    }

    for (int32_t i = 0; i < n; i++)
    {
        if (ibox_is_char(child[i].cibox))
        {
            if (w.back().empty() || !eis_str(w.back().back().get()))
                w.back().emplace_back(emake_str("", 0));
            std::string s = estr_mkstring(w.back().back().get());
            stdstring_pushback_char16(s, ibox_type(child[i].cibox)&65535);
            w.back().back().reset(emake_str(s));
        }
        else
        {
            boxbase * c = ibox_to_ptr(child[i].cibox);
            if (c->get_type() == BNTYPE_ROWBEGIN)
            {
                w.push_back(std::vector<wex>());
            }
            else if (c->get_type() == BNTYPE_ROWEND)
            {
                assert(!w.empty());
                ex e = finish_row(w.back());
                w.pop_back();
                w.back().emplace_back(e);
            }
            else
            {
                w.back().emplace_back(c->get_ex());
            }
        }
    }

    assert(w.size() == 1);
    return finish_row(w.back());
}

void rowbegin::get_cursor(aftransform * T)
{
    assert(false);
}

void rowend::get_cursor(aftransform * T)
{
    assert(false);
}

void rowbox::get_cursor(aftransform * T)
{
    assert(cursor_a < child.size());
    if (cursor_b >= child.size())
    {
        ibox_to_ptr(child[cursor_a].cibox)->get_cursor(T);
    }
    else
    {
        T->orig_x = 0;
        T->orig_y = ibox_centery(child[cursor_a].cibox);
        T->theta = 0;
        T->cos_theta = 1.0;
        T->sin_theta = 0.0;
    }
    T->orig_x += child[cursor_a].offx;
    T->orig_y += child[cursor_a].offy;
}


void rowbegin::measure(boxmeasurearg ma)
{
    assert(false);
}

void rowend::measure(boxmeasurearg ma)
{
    assert(false);
}

/*
struct rowoptions {
    boxmeasurearg ma;
    double default_sizex, default_sizey, default_centery;
    int32_t begin_idx;

    rowoptions(const boxmeasurearg & other, double default_sizex_, double default_sizey_, double default_centery_, int32_t begin_idx_)
        : ma(other),
          default_sizex(default_sizex_),
          default_sizey(default_sizey_),
          default_centery(default_centery_),
          begin_idx(begin_idx_)
    {}
};
*/

void rowbox::measure(boxmeasurearg mma)
{
//std::cout << "printing box: " << std::endl;
//    print(0,0,0);
//printf("rowbox  font: %08lx\n", mma.style->get_font());
//printf("   fontcolor: %08lx\n", mma.style->get_fontcolor());

    if (flags & BNFLAG_MEASURED)
        return;

//    mma.style->push_options(style.get());

    double maxabove, maxbelow, accwidth, accheight, maxwidth, maxheight;
//    rasterfont * font = fontint_to_fontp(fi);
//    int default_sizey = font->default_sizey();
//    int default_centery = font->default_centery();
//    int default_sizex = font->default_sizex();

//printf("fi: %d\n",fi);
    fs = mma.style->get_font();
    fcolor = mma.style->get_fontcolor();
//    double default_sizey = fontsize_default_sizey(fs);
//    double default_centery = fontsize_default_centery(fs);
//    double default_sizex = fontsize_default_sizex(fs);


//printf("rowbox measure fs: %08lx\n", fs);


//printf("default_sizex %f\n",default_sizex);
//printf("default_sizey %f\n",default_sizey);
//printf("default_centery %f\n",default_centery);

    mma.wrap_width = std::max(mma.wrap_width, 4*mma.style->get_default_sizex());
    mma.level++;

    flags |= BNFLAG_MEASURED;
    flags &= ~BNFLAG_COLORED;

//printf("boxnode_measure row\n"); boxnode_print(nullptr, box(us), 0);



    std::vector<int32_t> ostack;
    std::vector<double> maxabove_between;
    std::vector<double> maxbelow_between;
    blexer L(child.size() + 1);

//printf("before 1st scan:\n"); boxnode_print(NULL, Us, 0);


    /* first scan - find lexical types, position brackets, and position supers */
    ostack.push_back(-1);
    maxabove_between.push_back(mma.style->get_default_centery());
    maxbelow_between.push_back(mma.style->get_default_sizey() - mma.style->get_default_centery());

    for (int32_t i = 0; i < child.size(); i++)
    {
        assert(ostack.size() > 0);
        assert(ostack.size() == maxabove_between.size());
        assert(ostack.size() == maxbelow_between.size());

        if (ibox_is_char(child[i].cibox))
        {
            uint32_t a = ibox_type(child[i].cibox)&65535;
            uint32_t b;
            int32_t j = ostack.back();

            double a_fsizex = fontsize_char_sizex(mma.style->get_font(), a);
            double a_fsizey = fontsize_char_sizey(mma.style->get_font(), a);
            double a_fcentery = fontsize_char_centery(mma.style->get_font(), a);

//std::cout << "character a: " << char(a) << "  sizes: " << a_fsizex << ", " << a_fsizey << ", " << a_fcentery << std::endl;

            maxabove_between.back() = std::max(maxabove_between.back(), a_fcentery);
            maxbelow_between.back() = std::max(maxbelow_between.back(), a_fsizey - a_fcentery);

            /* the guess returned by add_char is the correct lextype for brackets */
            a += 65536*L.add_char(a);

            if (a == '(' + 65536*lextype_parenth_open ||
                a == '{'+ 65536*lextype_curly_open ||
                a == '['+ 65536*lextype_bracket_open ||
                a == CHAR_LeftDoubleBracket + 65536*lextype_bracket_open)
            {
                ostack.push_back(i);
                maxabove_between.push_back(a_fcentery);
                maxbelow_between.push_back(a_fsizey - a_fcentery);
            }
            else if (j >= 0 && ((b = '(', a == ')' + 65536*lextype_parenth_close) ||
                                (b = '{', a == '}' + 65536*lextype_curly_close) ||
                                (b = '[', a == ']' + 65536*lextype_bracket_close) ||
                                (b = CHAR_LeftDoubleBracket, a == CHAR_RightDoubleBracket + 65536*lextype_bracket_close))
                            && (ibox_is_char(child[j].cibox, b)))
            {
                // we have found a matching bracket b at index j
                if (maxabove_between.back() > mma.style->get_default_centery() ||
                    maxbelow_between.back() > mma.style->get_default_sizey() - mma.style->get_default_centery())
                {
                    double half = std::max(maxabove_between.back(), maxbelow_between.back());
                    if (a == '}' + 65536*lextype_curly_close ||
                        a == ']' + 65536*lextype_bracket_close ||
                        a == CHAR_RightDoubleBracket + 65536*lextype_bracket_close)
                    {
                        half = 0.375*half + 0.625*mma.style->get_default_centery();
                    }
                    else
                    {
                        half = 0.75*half + 0.25*mma.style->get_default_centery();
                    }

                    if (2*half > mma.style->get_default_sizey())
                    {
                        a_fcentery = half;
                        a_fsizey = 2*half;
                        child[j].cibox = iboxchar_set_sizes(child[j].cibox, a_fsizex, a_fsizey, a_fcentery);
                    }
                }
                assert(maxabove_between.size() == maxbelow_between.size());
                assert(maxabove_between.size() > 1);
                ostack.pop_back();
                double above = maxabove_between.back(); maxabove_between.pop_back();
                double below = maxbelow_between.back(); maxbelow_between.pop_back();
                maxabove_between.back() = std::max(maxabove_between.back(), above);
                maxbelow_between.back() = std::max(maxbelow_between.back(), below);
            }
            // set size of a:
            child[i].cibox = iboxchar_set_sizes(child[i].cibox, a_fsizex, a_fsizey, a_fcentery);
        }
        else if (ibox_to_ptr(child[i].cibox)->get_type() == BNTYPE_NULLER)
        {
            L.add_newline();
        }
        else if (ibox_to_ptr(child[i].cibox)->get_type() == BNTYPE_ROWBEGIN)
        {
            rowbegin * cc = dynamic_cast<rowbegin*>(ibox_to_ptr(child[i].cibox));
            cc->sizex = 0;
            cc->sizey = mma.style->get_default_sizey();
            cc->centery = mma.style->get_default_centery();
//            mma.style->push_options(cc->style.get());
            L.add_box(BNTYPE_ROWBEGIN);
        }
        else if (ibox_to_ptr(child[i].cibox)->get_type() == BNTYPE_ROWEND)
        {
//            mma.style->pop_options();
            L.add_box(BNTYPE_ROWEND);
        }
        else
        {
            boxbase* c = ibox_to_ptr(child[i].cibox);
            c->measure(mma);

            L.add_box(c->get_type());

            if (c->get_type() == BNTYPE_SUPER)
            {
                superscriptbox* s = dynamic_cast<superscriptbox*>(c);
                double prev_centery = i > 0 ? ibox_centery(child[i - 1].cibox) : mma.style->get_default_centery();
                double prev_sizey   = i > 0 ? ibox_sizey(child[i - 1].cibox) : mma.style->get_default_sizey();
                s->super.offx = 0;
                s->sizex = s->super.cbox->sizex;
                s->super.offy = 0;
                double cy1 = std::max(prev_centery, s->super.cbox->sizey - 1.0/16*(prev_sizey - prev_centery));
                double cy2 = prev_centery + s->super.cbox->centery;
                double alpha = prev_sizey/mma.style->get_default_sizey();
                alpha = std::max(alpha - 1, double(0));
                alpha = alpha/(alpha + 1);
                s->centery = (1-alpha)*cy1 + (alpha)*cy2;
                s->sizey = s->centery + (prev_sizey - prev_centery);
            }
            else if (c->get_type() == BNTYPE_SUB)
            {
                subscriptbox* s = dynamic_cast<subscriptbox*>(c);
                double prev_centery = i > 0 ? ibox_centery(child[i - 1].cibox) : mma.style->get_default_centery();
                double prev_sizey   = i > 0 ? ibox_sizey(child[i - 1].cibox) : mma.style->get_default_sizey();
                s->sub.offx = 0;
                s->sub.offy = std::max(prev_sizey - s->sub.cbox->sizey, 0.75*prev_centery);
                s->sizex = s->sub.cbox->sizex;
                s->sizey = s->sub.offy + s->sub.cbox->sizey;
                s->centery = prev_centery;
            }
            else if (c->get_type() == BNTYPE_SUBSUPER)
            {
                subsuperscriptbox* s = dynamic_cast<subsuperscriptbox*>(c);
                double prev_centery = mma.style->get_default_centery();
                double prev_sizey = mma.style->get_default_sizey();
                int32_t slant = 0;
                if (i > 0)
                {
                    prev_sizey   = ibox_sizey(child[i - 1].cibox);
                    prev_centery = ibox_centery(child[i - 1].cibox);
                    int32_t prev_type = ibox_type(child[i - 1].cibox);
                    /* adjust limits of integral */
                    if (prev_type >= 0 && (prev_type&65535) == CHAR_Integral)
                    {
                        slant = fontsize_size(fs)/3;
                    }
                }
                s->super.offx = slant;
                s->sub.offx = -slant;
                s->super.offy = std::min(-prev_centery, -s->super.cbox->sizey);
                s->sub.offy = std::max(prev_sizey - prev_centery - s->sub.cbox->sizey, 0.0);
                s->sizex = std::max(s->sub.cbox->sizex, s->super.cbox->sizex);
                s->sizey = s->sub.offy + s->sub.cbox->sizey - s->super.offy;
                s->centery = -s->super.offy;
                s->sub.offy -= s->super.offy;
                s->super.offy = 0;
            }

            maxabove_between.back() = std::max(maxabove_between.back(), c->centery);
            maxbelow_between.back() = std::max(maxbelow_between.back(), c->sizey - c->centery);
        }
    }

    /* second scan - place extra horiz space between appropriate lexical types */
    if (mma.style->get_autospacing())
    {
        int32_t us_lextype, prev_lextype = lextype_unknown;
        int32_t us_type, prev_type = 0;
        int32_t prev_is_expr = 0;
        for (int32_t i = 0; i < child.size(); i++)
        {
            us_lextype = L.type[i];
            if (ibox_is_char(child[i].cibox))
            {
                if (ibox_is_imm(child[i].cibox))
                {
                    us_type = iboximm_type(child[i].cibox)&65535;
                    child[i].cibox = iboximm_addlextype(child[i].cibox, us_lextype);
                }
                else
                {
                    us_type = ibox_to_ptr(child[i].cibox)->get_type()&65535;
                    ibox_to_ptr(child[i].cibox)->set_type(us_type + 65536*us_lextype);
                }
            }
            else
            {
                us_type = ibox_to_ptr(child[i].cibox)->get_type();            
            }

            us_lextype = us_lextype & 255;

            int32_t e = ((prev_lextype == lextype_comma) && (us_lextype != lextype_whitespace)) ? 4 : 0;

            switch (us_lextype)
            {
                case lextype_unknown:
                    break;

                case lextype_number:
                case lextype_symbol:
                    break;


                case lextype_curly_open:
                case lextype_parenth_open:
                    e = std::max(e, prev_is_expr);
                    prev_is_expr = 0;
                    break;

                case lextype_bracket_open:
                    prev_is_expr = 0;
                    break;

                case lextype_expr:
                    if (prev_is_expr && us_type != BNTYPE_SUB &&
                                        us_type != BNTYPE_SUBSUPER &&
                                        us_type != BNTYPE_SUPER)
                    {
                        e = 4;
                    }
                    prev_is_expr = 4;
                    break;

                case lextype_slot_1st:
                case lextype_blank_1st:
                case lextype_string_1st:
                case lextype_number_1st:
                case lextype_symbol_1st:
                case lextype_pattern_1st:
                    e = std::max(e, prev_is_expr);
                    prev_is_expr = 4;
                    break;

                case lextype_oppost_1st:
                    prev_is_expr = 4;
                    break;
                case lextype_oppre_1st:
                    e = std::max(e, prev_is_expr);
                    prev_is_expr = 0;
                    break;
                case lextype_opinf_1st:
                    prev_is_expr = (us_type == '/') ? 3 :
                                   (us_type == '*') ? 3 :
                                   (us_type == '^') ? 2 :
                                                      4;
                    if (us_type != ';')
                        e = std::max(e, prev_is_expr);
                    break;

                case lextype_opinf:
                    if (prev_is_expr < 4)
                    {
                        assert(i > 0);
                        L.extrax[i - 1] = 4;
                        prev_is_expr = 4;
                    }
                case lextype_oppost:
                case lextype_oppre:
                    break;

                case lextype_curly_close:
                case lextype_parenth_close:
                case lextype_bracket_close:
                    prev_is_expr = 4;
                    break;

                case lextype_comma:
                    prev_is_expr = 0;
                    break;

                case lextype_message_name:
                case lextype_string:
                case lextype_blank:
                case lextype_slot:
                case lextype_comment:
                case lextype_comment_1st:
                case lextype_whitespace:
                    break;
                default:
                    std::cout << "unknown lextype: " << us_lextype << std::endl;
                    assert(false);
            }

            if (prev_lextype == lextype_unknown || prev_lextype == lextype_whitespace)
            {
                e = 0;
            }

            L.extrax[i] = e;

            prev_lextype = us_lextype;
            prev_type = us_type;
        }
    }
    else
    {
        int32_t us_lextype;
        for (int32_t i = 0; i < child.size(); i++)
        {
            us_lextype = L.type[i];
            if (ibox_is_char(child[i].cibox))
            {
                if (ibox_is_imm(child[i].cibox))
                {
                    child[i].cibox = iboximm_addlextype(child[i].cibox, us_lextype);
                }
                else
                {
                    ibox_to_ptr(child[i].cibox)->set_type((ibox_to_ptr(child[i].cibox)->get_type() & 65535) + 65536*us_lextype);
                }
            } else {
                //us_type = bptr_type(child);
            }

            L.extrax[i] = 0;
        }
    }

//printf("NEW before 3rd scan:\n");

    /* third scan - find line breaks and place children into row */
    accheight = 0;
    double accheightsave = accheight;
//std::cout << "accheight: " << accheight << std::endl;
    maxwidth = 1;
    maxabove = mma.style->get_default_centery();
    maxbelow = mma.style->get_default_sizey() - mma.style->get_default_centery();
    int32_t a = 0, b = 0;
    int32_t checklinecount=0;
    double accwidth_restart = mma.style->get_autospacing() ? mma.style->get_default_sizex() : 0;

//printf("processing new row\n");
    do {
//printf("a = %d\n",a);
        checklinecount++;
        // find the delimiters [a,b) of this line
        while (b + 1 < child.size())
        {
            if (ibox_type(child[b].cibox) == BNTYPE_NULLER)
            {
                break;
            }
            b++;
        }

        // process chars in [a,b)
        int32_t child_type;
        double child_sizex, child_sizey, child_centery;
        int32_t c = a;
        int32_t d = a;
        accwidth = 0;
        bool processed_end = false;
        while (d < b)
        {//printf("looking at [%d, %d)\n",c,d); printf("accwidth %d\n",accwidth); SleepUS(50000);
            child_type = ibox_type(child[d].cibox);
            child_sizex = ibox_sizex(child[d].cibox);
            child_sizey = ibox_sizey(child[d].cibox);
            child_centery = ibox_centery(child[d].cibox);

            child[d].offx = accwidth;
            //us->array[d].offy = accheight + maxabove - child_centery;

            accwidth += child_sizex + L.extrax[d + 1]*0.25*mma.style->get_default_sizex();

            if (accwidth > mma.wrap_width && c < d)
            {
                int32_t j = d;
                int32_t smallest = L.penalty[j] + (d-j)*(d-j);
                for (int32_t i = d - 1; i > c; i--)
                {
                    if (smallest > L.penalty[i] + (d-i)*(d-i))
                    {
                        smallest = L.penalty[i] + (d-i)*(d-i);
                        j = i;
                    }
                }
                d = j;
                assert(c < d);
                maxabove = mma.style->get_default_centery();
                maxbelow = mma.style->get_default_sizey() - mma.style->get_default_centery();

                for (int32_t i = c; i < d; i++)
                {
                    child_sizey = ibox_sizey(child[i].cibox);
                    child_centery = ibox_centery(child[i].cibox);
                    maxabove = std::max(maxabove, child_centery);
                    maxbelow = std::max(maxbelow, child_sizey - child_centery);
                }
                for (int32_t i = c; i < d; i++)
                {
                    child[i].offy = accheight + maxabove - ibox_centery(child[i].cibox);
                }
                accheight += maxabove + maxbelow;
                accheightsave = accheight;
                accheight += 3.0/32*mma.style->get_default_sizey();

                maxwidth = std::max(maxwidth, child[d].offx);

                processed_end = true;
                if (d >= b)
                {
                    // process ending
                    assert(d == b);
                    child[d].offx = accwidth;
                    child[d].offy = accheight + maxabove - ibox_centery(child[d].cibox);
                }

                accwidth = accwidth_restart;
                c = d;
            }
            else
            {
                processed_end = false;
                d++;
            }
        }
        assert(d == b);
        if (!processed_end)
        {
                maxabove = mma.style->get_default_centery();
                maxbelow = mma.style->get_default_sizey() - mma.style->get_default_centery();

                for (int32_t i = c; i < d; i++)
                {
                    child_sizey = ibox_sizey(child[i].cibox);
                    child_centery = ibox_centery(child[i].cibox);
                    maxabove = std::max(maxabove, child_centery);
                    maxbelow = std::max(maxbelow, child_sizey - child_centery);
                }
                for (int32_t i = c; i < d; i++)
                {
                    child[i].offy = accheight + maxabove - ibox_centery(child[i].cibox);
                }

                // process ending
                child[d].offx = accwidth;
                child[d].offy = accheight + maxabove - ibox_centery(child[d].cibox);

                accheight += maxabove + maxbelow;
                accheightsave = accheight;
                accheight += 3.0/32*mma.style->get_default_sizey();
                maxwidth = std::max(maxwidth, child[d].offx);
        }

        // goto next line
        a = ++b;
        if (a < child.size())
        {
            accheightsave = accheight;
            accheight += 5.0/32 * mma.style->get_default_sizey();
        }

    } while (a < child.size());
    // process ending
//        child = us->array[us->len-1].child;
//        us->array[us->len-1].offx = accwidth;
//        us->array[us->len-1].offy = accheight + maxabove - child->centery;

//printf("maxwidth: %d  accwidth: %d  maxabove: %d\n",maxwidth,accwidth,maxabove);
//printf("checklinecount: %d\n",checklinecount);

    sizex = maxwidth;
    sizey = accheightsave;
    centery = child[0].offy + ibox_centery(child[0].cibox);


//    centery = checklinecount == 1 ? maxabove : 0.5*accheight;


//    mma.style->pop_options();

//printf("measure row 3 done\n");
}



/* draw ******************************/

void rowbegin::draw_pre(boxdrawarg da)
{
    assert(false);
}

void rowend::draw_pre(boxdrawarg da)
{
    assert(false);
}

void rowbox::draw_pre(boxdrawarg da)
{
//std::cout << "rowbox::draw_pre " << da.tostring() << std::endl;

    if (cursor_b >= child.size())
    {
        assert(cursor_a < child.size());
        ibox_to_ptr(child[cursor_a].cibox)->draw_pre(boxdrawarg(da, child[cursor_a].offx, child[cursor_a].offy));
        return;
    }
    else if (cursor_a == cursor_b)
    {
        int32_t idx  = cursor_a;
        std::vector<int32_t> ostack;
        std::vector<int32_t> commapos;
        int32_t match = 0;
        int32_t right_pi = -1;
        int32_t right_x = 0;
        ostack.clear();
        for (int32_t count = 5000; count >= 0; count--)
        {
            int32_t t = ibox_type(child[idx].cibox);
            if (t >= 0)
            {
                int32_t y, x = t & 0x0FFFFFF;
                if (x == ',' + 65536*lextype_comma && ostack.empty())
                {
                    commapos.push_back(idx);
                }
                else if (   (y =  '(' + 65536*lextype_parenth_open,
                             x == ')' + 65536*lextype_parenth_close)
                         || (y =  '{' + 65536*lextype_curly_open,
                             x == '}' + 65536*lextype_curly_close)
                         || (y =  '[' + 65536*lextype_bracket_open,
                             x == ']' + 65536*lextype_bracket_close)
                         || (y =  CHAR_LeftDoubleBracket + 65536*lextype_bracket_open,
                             x == CHAR_RightDoubleBracket + 65536*lextype_bracket_close))
                {
                    if (ostack.empty())
                    {
                        right_pi =idx;
                        right_x = x;
                        break;
                    }
                    if (y == ostack.back())
                    {
                        ostack.pop_back();
                    }
                }
                else if (   x == '(' + 65536*lextype_parenth_open
                         || x == '{' + 65536*lextype_curly_open
                         || x == '[' + 65536*lextype_bracket_open
                         || x == CHAR_LeftDoubleBracket + 65536*lextype_bracket_open)
                {
                    ostack.push_back(x);
                }
            }
            ++idx;
            if (idx >= child.size())
            {
                break;
            }
        }

        idx = cursor_a;
        match = 0;
        int32_t left_pi = -1;
        int32_t left_x = 0;
        ostack.clear();
        for (int32_t count = 5000; count >= 0; count--)
        {
            --idx;
            if (idx < 0)
            {
                break;
            }
            int32_t t = ibox_type(child[idx].cibox);
            if (t >= 0)
            {
                int32_t y, x = t & 0x0FFFFFF;
                if (x == ',' + 65536*lextype_comma && ostack.empty())
                {
                    commapos.push_back(idx);
                }
                else if (   (y =  ')' + 65536*lextype_parenth_close,
                             x == '(' + 65536*lextype_parenth_open)
                         || (y =  '}' + 65536*lextype_curly_close,
                             x == '{' + 65536*lextype_curly_open)
                         || (y =  ']' + 65536*lextype_bracket_close,
                             x == '[' + 65536*lextype_bracket_open)
                         || (y =  CHAR_RightDoubleBracket + 65536*lextype_bracket_close,
                             x == CHAR_LeftDoubleBracket + 65536*lextype_bracket_open))
                {
                    if (ostack.size() < 1)
                    {
                        left_pi = idx;
                        left_x = x;
                        break;
                    }
                    if (y == ostack.back())
                    {
                        ostack.pop_back();
                    }
                }
                else if (   x == ')' + 65536*lextype_parenth_close
                         || x == '}' + 65536*lextype_curly_close
                         || x == ']' + 65536*lextype_bracket_close
                         || x == CHAR_RightDoubleBracket + 65536*lextype_bracket_close)
                {
                    ostack.push_back(x);
                }
            }
        }

        uint32_t bracketcolor = da.nb->cBracketMismatch;
        if (   (    left_x == '[' + 65536*lextype_bracket_open
                && right_x == ']' + 65536*lextype_bracket_close)
            || (    left_x == '(' + 65536*lextype_parenth_open
                && right_x == ')' + 65536*lextype_parenth_close)
            || (    left_x == '{' + 65536*lextype_curly_open
                && right_x == '}' + 65536*lextype_curly_close)
            || (    left_x == CHAR_LeftDoubleBracket + 65536*lextype_bracket_open
                && right_x == CHAR_RightDoubleBracket + 65536*lextype_bracket_close))
        {
            bracketcolor = da.nb->cBracketMatch;
            for (auto& i : commapos)
            {
                double sib_sizex = ibox_sizex(child[i].cibox);
                double sib_sizey = ibox_sizey(child[i].cibox);
                drawtrect(da.globx + child[i].offx, da.globx + child[i].offx + sib_sizex,
                          da.globy + child[i].offy + sib_sizey/2, da.globy + child[i].offy + sib_sizey, bracketcolor, da.T);
            }
        }
        if (left_x != 0)
        {
            double sib_sizex = ibox_sizex(child[left_pi].cibox);
            double sib_sizey = ibox_sizey(child[left_pi].cibox);
            drawtrect(da.globx + child[left_pi].offx, da.globx + child[left_pi].offx + sib_sizex,
                      da.globy + child[left_pi].offy, da.globy + child[left_pi].offy + sib_sizey, bracketcolor, da.T);
        }
        if (right_x != 0)
        {
            double sib_sizex = ibox_sizex(child[right_pi].cibox);
            double sib_sizey = ibox_sizey(child[right_pi].cibox);
            drawtrect(da.globx + child[right_pi].offx, da.globx + child[right_pi].offx + sib_sizex,
                      da.globy + child[right_pi].offy, da.globy + child[right_pi].offy + sib_sizey, bracketcolor, da.T);
        }
    }
    else
    {
        double selstart = std::min(cursor_a, cursor_b);
        double selend = std::max(cursor_a, cursor_b);
        for (int32_t i = selstart; i < selend; i++)
        {
            double child_sizex = ibox_sizex(child[i].cibox);
            double child_sizey = ibox_sizey(child[i].cibox);
            double usx = da.globx + child[i].offx;
            double usy = da.globy + child[i].offy;
            if (child_sizex > 0)
            {
                drawtrect(usx, usx + child_sizex,
                          usy, usy + child_sizey, da.nb->cSelectionBackground, da.T);
            }
        }
    }

    // draw the cursor
    double default_sizey = fontsize_default_sizey(fs);
    double default_centery = fontsize_default_centery(fs);
    int32_t us_type = ibox_type(child[cursor_a].cibox);
    double us_sizex = ibox_sizex(child[cursor_a].cibox);
    double us_sizey = ibox_sizey(child[cursor_a].cibox);
    double us_centery = ibox_centery(child[cursor_a].cibox);
    double cursor_y;
    double cursor_sizey;
    if (cursor_a > 0 && ibox_type(child[cursor_a - 1].cibox) != BNTYPE_NULLER)
    {
        double prev_sizey, prev_centery;
        double us_offy = child[cursor_a].offy;
        double prev_offy = child[cursor_a - 1].offy;
        prev_sizey = ibox_sizey(child[cursor_a - 1].cibox);
        prev_centery = ibox_centery(child[cursor_a - 1].cibox);
        if (us_type != BNTYPE_NULLER)
        {
            cursor_y = (us_centery - prev_centery)/2;
            cursor_sizey = (prev_sizey + us_sizey)/2;
        }
        else
        {
            cursor_y = us_centery - prev_centery;
            cursor_sizey = prev_sizey;
        }
    }
    else
    {
        if (us_type != BNTYPE_NULLER)
        {
            cursor_y = 0;
            cursor_sizey = us_sizey;
        }
        else
        {
            cursor_y = us_centery - default_centery;
            cursor_sizey = default_sizey;
        }
    }

    double usx = da.globx + child[cursor_a].offx;
    double usy = da.globy + child[cursor_a].offy;
    double f = fontsize_size(fs);
    if (us_type == BNTYPE_ROT)
    {
        double Ax = usx;
        double Ay = usy;
        double Bx = usx + us_sizex;
        double By = usy;
        double Cx = usx;
        double Cy = usy + us_sizey;
        double Dx = usx + us_sizex;
        double Dy = usy + us_sizey;
        drawtline(Ax, Ay, Bx, By, f*0.0625, fcolor, da.T);
        drawtline(Bx, By, Dx, Dy, f*0.0625, fcolor, da.T);
        drawtline(Dx, Dy, Cx, Cy, f*0.0625, fcolor, da.T);
        drawtline(Cx, Cy, Ax, Ay, f*0.125, da.nb->cCursor, da.T);
    }
    else
    {
        drawtline(usx, usy + cursor_y,
                  usx, usy + cursor_y + cursor_sizey,
                  f*0.25, da.nb->cCursor, da.T);
    }

    return;
}


void rowbegin::draw_main(boxdrawarg da)
{
    assert(false);
}

void rowend::draw_main(boxdrawarg da)
{
    assert(false);
}

void rowbox::draw_main(boxdrawarg da)
{
//std::cout << "rowbox::draw_main " << da.tostring() << std::endl;

//    double default_sizey = fontsize_default_sizey(fs);
//    double default_centery = fontsize_default_centery(fs);

    // whole box is off screen TODO: rotation
    if (   (da.globx >= glb_image.sizex || da.globx + sizex <= 0)
        && (da.globy >= glb_image.sizey || da.globy + sizey <= 0) )
    {
        return;
    }

    assert(!child.empty());

    // binary search to find last y coor y1 <= 0
    int32_t y1, y1i;
    {
        int32_t low, mid, hi;
        int32_t lowy, midy, hiy;

        low = 0;
        hi = child.size() - 1;
        lowy = da.globy + child[low].offy + ibox_centery(child[low].cibox);
        hiy = da.globy + child[hi].offy + ibox_centery(child[hi].cibox);

        if (lowy >= 0) {
            y1i = low;
            y1 = lowy;
        } else if (hiy <= 0) {
            y1i = hi;
            y1 = hiy;
        } else {
            // lowy < 0 and hiy > 0
            while (low < hi) {
                mid = (low + hi)/2;
                midy = da.globy + child[mid].offy + ibox_centery(child[mid].cibox);
                if (midy <= 0) {
                    if (mid == low)
                        break;
                    low = mid;
                    lowy = midy;
                    if (midy == 0)
                        break;
                } else {
                    hi = mid;
                    hiy = midy;
                }
            }
            y1i = low;
            y1 = lowy;
        }
    }

    // binary search to find last y coor y0 < y1;
    int32_t y0, y0i;
    {
        int32_t low, mid, hi;
        int32_t lowy, midy, hiy;

        low = 0;
        hi = y1i;
        lowy = da.globy + child[low].offy + ibox_centery(child[low].cibox);
        hiy = da.globy + child[hi].offy + ibox_centery(child[hi].cibox);

        if (lowy >= y1) {
            y0i = low;
            y0 = lowy;
        } else if (hiy < y1) {
            y0i = hi;
            y0 = hiy;
        } else {
            // lowy < y1 and hiy > y1
            while (low < hi) {
                mid = (low + hi)/2;
                midy = da.globy + child[mid].offy + ibox_centery(child[mid].cibox);
                if (midy < y1) {
                    if (mid == low)
                        break;
                    low = mid;
                    lowy = midy;
                } else {
                    hi = mid;
                    hiy = midy;
                }
            }
            y0i = low;
            y0 = lowy;
        }
    }

    // set start and stop indices for drawing
    int32_t starti = y0i, stopi;
    int32_t lasty = 0;
    for (stopi = starti; stopi < child.size(); stopi++)
    {
        if (lasty > int32_t(glb_image.sizey) && da.globy + child[stopi].offy > lasty)
        {
            break;
        }
        if (lasty <= int32_t(glb_image.sizey))
        {
            lasty = da.globy + child[stopi].offy + ibox_centery(child[stopi].cibox);
        }
    }

    int32_t selstart = 0;
    int32_t selend = 0;

    if (da.dflags & DFLAG_SELECTION)
    {
        selend = child.size();
    }
    else if (((da.dflags & DFLAG_IGNORESEL) == 0) && cursor_b < child.size())
    {
        selstart = std::min(cursor_a, cursor_b);
        selend = std::max(cursor_a, cursor_b);
    }

    int32_t asdf = cursor_b >= child.size() ? cursor_a : -1;

    // lets draw
    if (da.dflags & DFLAG_SCOLOR)
    {
        if (!(flags & BNFLAG_COLORED))
        {
            bsymer B;
            B.add_scolor(child.data(), 0, child.size() - 1);
            flags |= BNFLAG_COLORED;
        }
    }

    std::vector<std::pair<uint32_t, uint32_t>> font_stack;

    font_stack.push_back(std::pair<uint32_t, uint32_t>(fs, fcolor));
    std::pair<uint32_t, uint32_t> * cur_font = &font_stack.back();

    int32_t end_balance = 0;
    for (int32_t i = 0; i < starti; i++)
    {
        if (ibox_is_char(child[i].cibox))
            continue;

        boxbase * cb = ibox_to_ptr(child[i].cibox);

        if (cb->get_type() == BNTYPE_ROWBEGIN)
        {
            rowbegin * cbb = dynamic_cast<rowbegin*>(cb);
            font_stack.push_back(std::pair<uint32_t, uint32_t>(cbb->fs, cbb->fcolor));
            cur_font = &font_stack.back();
        }
        else if (cb->get_type() == BNTYPE_ROWEND)
        {
            font_stack.pop_back();
            cur_font = &font_stack.back();
        }
    }

    int32_t contx = 0;
    int32_t conty = sizey;
    for (int32_t i = starti; i < stopi; i++)
    {
        iboxarrayelem * c = child.data() + i;
        if (ibox_is_char(c->cibox))
        {
            int32_t child_type = ibox_type(c->cibox);
            double child_sizex = ibox_sizex(c->cibox);
            double child_sizey = ibox_sizey(c->cibox);
            double child_centery = ibox_centery(c->cibox);
            uint32_t color;

            if (selstart <= i && i < selend)
            {
                color = da.nb->cSelectionForeground;
            }
            else if (da.dflags & DFLAG_SCOLOR)
            {
                if (((child_type>>16)&255) == lextype_symbol_1st ||
                    ((child_type>>16)&255) == lextype_pattern_1st ||
                    ((child_type>>16)&255) == lextype_symbol)
                {
                    color = da.nb->pallet1[lextype_MAX + ((child_type>>24)&255)];
                }
                else
                {
                    color = da.nb->pallet1[(child_type>>16)&255];
                }
            }
            else
            {
                color = cur_font->second;
            }

            drawtchar(cur_font->first, child_type&65535, child_sizex, child_sizey, da.globx + c->offx, da.globy + c->offy, color, da.T);
            int32_t newy = c->offy + child_centery;
//                if ((da.dflags & DFLAG_NLINE) && newy > conty + 2){bitmap_draw_rding(dingbats.continuation, image, globx + contx, globy + conty - dingbats.continuation.centery, dingbats.continuation.sizey, pallet1[lextype_unknown]);}
            contx = c->offx + child_sizex;
            conty = newy;
        }
        else if (ibox_to_ptr(child[i].cibox)->get_type() == BNTYPE_ROWBEGIN)
        {
            rowbegin * cbb = dynamic_cast<rowbegin*>(ibox_to_ptr(child[i].cibox));
            font_stack.push_back(std::pair<uint32_t, uint32_t>(cbb->fs, cbb->fcolor));
            cur_font = &font_stack.back();
        }
        else if (ibox_to_ptr(child[i].cibox)->get_type() == BNTYPE_ROWEND)
        {
            font_stack.pop_back();
            cur_font = &font_stack.back();
        }
        else
        {
            boxbase * cboxp = ibox_to_ptr(c->cibox);
            if (cboxp->get_type() == BNTYPE_NULLER)
            {
                contx = 0;
                conty = sizey;
            }
            else
            {
                cboxp->draw_main(boxdrawarg(da, c->offx, c->offy, (selstart <= i && i < selend) ? DFLAG_SELECTION : i == asdf ? 0 : DFLAG_IGNORESEL));
                double newy = c->offy + cboxp->centery;
//                    if ((da.dflags & DFLAG_NLINE) && newy > conty + 2){bitmap_draw_rding(dingbats.continuation, image, globx + contx, globy + conty - dingbats.continuation.centery, dingbats.continuation.sizey, pallet1[lextype_unknown]);}
                contx = c->offx + cboxp->sizex;
                conty = newy;
            }
        }
    }
}

void rowbegin::draw_post(boxdrawarg da)
{
    assert(false);
}

void rowend::draw_post(boxdrawarg da)
{
    assert(false);
}

void rowbox::draw_post(boxdrawarg da)
{
//std::cout << "rowbox::draw_post " << da.tostring() << std::endl;

    if (cursor_b >= child.size())
    {
        assert(cursor_a < child.size());
        ibox_to_ptr(child[cursor_a].cibox)->draw_post(boxdrawarg(da, child[cursor_a].offx, child[cursor_a].offy));
        return;
    }
    else if (cursor_a == cursor_b)
    {
        double default_sizey = fontsize_default_sizey(fs);
        double default_centery = fontsize_default_centery(fs);
        int32_t us_type = ibox_type(child[cursor_a].cibox);

        if (da.dflags && DFLAG_SCOLOR)
        {
            double pk = cursor_a;
            double max_y = 0;
            if (!(ibox_is_char(child[pk].cibox) >= 0
                  && isletterchar(char16_t(ibox_type(child[pk].cibox)))))
            {
                while (pk > 0 && ibox_is_char(child[pk - 1].cibox)
                              && isletterchar(char16_t(ibox_type(child[pk - 1].cibox))))
                {
                    max_y = std::max(max_y, child[pk - 1].offy + ibox_sizey(child[pk - 1].cibox));
                    pk--;
                }
            }

            if (pk < cursor_a - 1)
            {
                std::string s;
                for (int32_t k = pk; k < cursor_a; k++)
                {
                    stdstring_pushback_char16(s, ibox_type(child[k].cibox)&65535);
                }
//    std::cout << "prefix " << s << std::endl;                    

                // Prefix search
                auto set_range = gs.char_set.equal_prefix_range(s);

                double oox, ox, oy = max_y;

                for(auto it = set_range.first; it != set_range.second; ++it)
                {
                    oox = child[pk].offx;
                    it.key(s);
                    double max_above = 0, max_below = 0;
                    ox = oox;
                    for (size_t si = 0; si < s.size();)
                    {
                        char16_t c;
                        si += readonechar16(c, s.data() + si);
                        double c_fsizex = fontsize_char_sizex(fs, c);
                        double c_fsizey = fontsize_char_sizey(fs, c);
                        double c_fcentery = fontsize_char_centery(fs, c);
                        max_above = std::max(max_above, c_fcentery);
                        max_below = std::max(max_below, c_fsizey - c_fcentery);
                        ox += c_fsizex;
                    }

                    drawbtrect(0.75, da.globx + oox, da.globx + ox,
                                     da.globy + oy, da.globy + oy + max_above + max_below, 0x00f0f0f0, da.T);

                    ox = oox;
                    for (size_t si = 0; si < s.size();)
                    {
                        char16_t c;
                        si += readonechar16(c, s.data() + si);
                        double c_fsizex = fontsize_char_sizex(fs, c);
                        double c_fsizey = fontsize_char_sizey(fs, c);
                        double c_fcentery = fontsize_char_centery(fs, c);

                        drawbtchar(0.75, fs, c, c_fsizex, c_fsizey, da.globx + ox, da.globy + oy + max_above - c_fcentery, 0x00101010, da.T);
                        ox += c_fsizex;
                    }
                    oy += max_above + max_below;
                }
            }
        }

        return;
    }
    else
    {
        return;
    }
}
