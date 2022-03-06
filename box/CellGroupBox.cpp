#include <stack>
#include "boxes.h"
#include "ex_parse.h"
#include "notebook.h"
#include "box_convert.h"

/* edit ************************************/

void cellgroupbox::print(size_t indent, double offx, double offy)
{
    for (size_t i = 0; i < indent; i++)
        printf("    ");
    printf("(%f,%f) cellgroup: size (%f,%f:%f)\n", offx, offy, sizex, sizey, centery);

    for (auto& c: child)
    {
        c.cbox->print(indent + 1, c.offx, c.offy);
    }
}

visitRet cellgroupbox::visit(visitArg m)
{
    switch (m)
    {
        case visitarg_InvalidateAll:
        {
//            flags &= ~(BNFLAG_MEASURED | BNFLAG_COLORED);
            for (auto i = child.begin(); i != child.end(); ++i)
            {
                i->cbox->visit(m);
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

boxbase * cellgroupbox::copy()
{
    size_t n = child.size();
    cellgroupbox * r = new cellgroupbox(n);
    for (size_t i = 0; i < n; i++)
    {
        r->child[i].cbox = child[i].cbox->copy();
    }
    return r;
}

ex cellgroupbox::get_ex()
{
    uex v; v.init_push_backr(gs.sym_sList.get(), child.size());
    for (int32_t i = 0; i < child.size(); i++)
        v.push_back(child[i].cbox->get_ex());
    return emake_node(gs.sym_sCellGroupData.copy(), v.release());
}

/* measure **********************************/

void cellgroupbox::measure(boxmeasurearg ma)
{
    double accheight = 0;
    double maxwidth = 0;
    for (auto & i : child)
    {
        i.cbox->measure(ma);
        i.offx = 0;
        i.offy = accheight;
        maxwidth = std::max(maxwidth, i.cbox->sizex);
        accheight += i.cbox->sizey;
    }

    sizex = maxwidth;
    sizey = accheight;
    centery = accheight/2;
}

/* draw *************************************/

void cellgroupbox::draw_pre(boxdrawarg da)
{
    return;
}

void cellgroupbox::draw_main(boxdrawarg da)
{
    _draw_cellgroup_bracket(this, da, false);
    for (auto & i : child)
    {
        i.cbox->draw_main(boxdrawarg(da, i.offx, i.offy));
    }
    return;
}

void cellgroupbox::draw_post(boxdrawarg da)
{
    return;
}
