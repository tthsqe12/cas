#include <stack>
#include "boxes.h"
#include "ex_parse.h"
#include "notebook.h"
#include "box_convert.h"

void cellbox::print(size_t indent, double offx, double offy)
{
    for (size_t i = 0; i < indent; i++)
        printf("    ");
    printf("(%f,%f) cell: size(%f,%f:%f)  measured: %d\n", offx, offy, sizex, sizey, centery, flags & BNFLAG_MEASURED);

    body.cbox->print(indent + 1, body.offx, body.offy);

    if (0 && label.cbox != nullptr)
        label.cbox->print(indent + 1, label.offx, label.offy);

    if (mexpr != nullptr)
        mexpr->print(indent + 1, 0, 0);
}


ex eremove_rule(er E, er L)
{
    uex e(ecopy(E));
    for (ulong i = 1; i <= elength(E); i++)
    {
        er ei = echild(E,i);
        if (!ehas_head_sym_length(ei, gs.sym_sRule.get(), 2))
            continue;
        if (!ex_same(echild(ei,1), L))
            continue;
        e.removechild(i);
        break;
    }
    return e.release();;
}



visitRet cellbox::visit(visitArg m)
{
    if (mexpr != nullptr)
        mexpr->visit(m);
    body.cbox->visit(m);
    if (label.cbox != nullptr)
        label.cbox->visit(m);
    return visitret_OK;
}

boxbase * cellbox::copy()
{
//    monobox* newmexpr = mexpr == nullptr ? nullptr : dynamic_cast<monobox*>(mexpr->copy());
    rowbox* newbody = dynamic_cast<rowbox*>(body.cbox->copy());
//    rowbox * newlabel = (label.cbox == nullptr) ? label.cbox : dynamic_cast<rowbox*>(label.cbox->copy());
//    rowbox* newlabel = nullptr;
    cellbox* r = new cellbox(newbody);
    r->options = options;
//    r->body.offx = body.offx;
//    r->body.offy = body.offy;
//    r->label.offx = label.offx;
//    r->label.offy = label.offy;
    return r;
}

void cellbox::key_copy(boxbase*&b)
{
    if (mexpr != nullptr)
    {
        mexpr->key_copy(b);
        return;
    }

    body.cbox->key_copy(b);
}

void cellbox::key_paste(boxbase*&b)
{
    if (mexpr != nullptr)
    {
        mexpr->key_paste(b);
        return;
    }

    body.cbox->key_paste(b);
}

void cellbox::insert_char(int32_t c)
{
    if (mexpr != nullptr)
    {
        mexpr->insert_char(c);
    }
    else
    {
        if (get_celllabelautodelete())
            options.remove(opt_CellLabel);

        body.cbox->insert_char(c);
    }
    flags &= ~BNFLAG_MEASURED;
    return;
}

moveRet cellbox::move(boxbase*&b, moveArg m)
{
    moveRet r;

    if (mexpr != nullptr)
        return mexpr->move(b, m);

    switch (m)
    {
        case movearg_Left:
        case movearg_ShiftLeft:
        case movearg_Right:
        case movearg_ShiftRight:
        case movearg_Up:
        case movearg_ShiftUp:
        case movearg_Down:
        case movearg_ShiftDown:
        case movearg_Tab:
        {
            r = body.cbox->move(b, m);
            assert(r == moveret_OK || r == moveret_End);
            return r;
        }
        case movearg_First:
        case movearg_Last:
        {
            r = body.cbox->move(b, m);
            assert(r == moveret_OK);
            return r;
        }
        case movearg_Home:
        case movearg_End:
        case movearg_CtrlSpace:
        {
            r = body.cbox->move(b, m);
            assert(r == moveret_OK || r == moveret_End);
            return r;
        }
        case movearg_Switch:
        {
            r = body.cbox->move(b, movearg_Switch);
            return moveret_OK;
        }
        default:
        {
            assert(false);
            return moveret_OK;
        }
    }
}

insertRet cellbox::insert(boxbase*&b, insertArg m)
{
    if (mexpr != nullptr)
        return mexpr->insert(b, m);

    if (get_celllabelautodelete())
        options.remove(opt_CellLabel);

    flags &= !BNFLAG_MEASURED;

    return body.cbox->insert(b, m);
}

removeRet cellbox::remove(boxbase*&b, removeArg m)
{
//std::cout << "cellbox::remove called" << std::endl;
    assert(b == nullptr);
    removeRet r;

    if (mexpr != nullptr)
    {
        r = mexpr->remove(b, m);
    }
    else
    {
        if (get_celllabelautodelete())
            options.remove(opt_CellLabel);

        r = body.cbox->remove(b, m);
    }

    assert(r == removeret_OK || r == removeret_End);

    flags &= !BNFLAG_MEASURED;
    return removeret_OK;
}

void cellbox::toggle_cell_expr()
{
    if (mexpr != nullptr)
    {
        install_live_psymbols LPS;
        eparser P;
        syntax_report sr;

        for (auto y = mexpr->array.begin(); y != mexpr->array.end(); ++y)
        {
            for (auto x = y->begin(); x != y->end(); ++x)
            {
                P.handle_rawchar(*x&65535);
                if (P.error)
                {
                    sr.handle_cstr_error(P.error, NULL, 0);
                    break;
                }
            }
            P.handle_newline();
            if (P.error)
            {
                sr.handle_row_error(P.error, NULL, 0);
                break;
            }
        }

        if (sr.error)
        {
            std::cerr << "syntax error" << std::endl;
            return;
        }
        P.handle_end();
        if (P.have_one_ex())
        {
            gs.send_psymbols();
            boxbase* newcell = boxbase_from_ex(P.estack[0].get());
            if (newcell->get_type() == BNTYPE_CELL)
            {
                cellbox* C = dynamic_cast<cellbox*>(newcell);
                std::swap(body.cbox, C->body.cbox);
                std::swap(options, C->options);
                flags &= ~BNFLAG_MEASURED;
                delete mexpr;
                mexpr = nullptr;
            }
            delete newcell;
            return;
        }
        else
        {
            std::cerr << "syntax error" << std::endl;
            return;
        }
    }
    else
    {
        wex e(get_ex());
        std::string s = ex_tostring_full(e.get());
        mexpr = new monobox(1, 0, 0, 0, 0);
        for (size_t i = 0; i < s.size(); )
        {
            char16_t c;
            i += readonechar16(c, s.c_str() + i);
            mexpr->array[0].push_back(c);
        }
    }
}

ex cellbox::get_ex()
{
    std::vector<wex> v;
    v.push_back(wex(body.cbox->get_ex()));
    options.get_ex(v);
    return emake_node(gs.sym_sCell.copy(), v);
}


void cellbox::get_cursor(aftransform * T)
{
    if (mexpr != nullptr)
    {
        mexpr->get_cursor(T);
        T->orig_x += bracket_offy;
        T->orig_y += bracket_offy;
    }
    else
    {
        body.cbox->get_cursor(T);
        T->orig_x += body.offx;
        T->orig_y += body.offy;
    }
}




void cellbox::measure_quick(stylestack * style)
{
    notebook * nb = style->parent_nb();

    options.set_restore_point();
    options.append_ex_list(style->get_opt(opt_CellOptions));


    for (ulong i = options.array.size(); i > 0 ; i--)
    {
        optType ot = options.array[i-1].option;
        if (ot == opt_CellLabel)
        {
            options.array[i-1].processed = true;
        }
        else if (ot == opt_CellLabelOptions)
        {
            options.array[i-1].processed = true;
        }
        else
        {
            options.array[i-1].processed = false;
        }
    }

    style->start_push_options();

    for (ulong i = options.array.size(); i > 0 ; i--)
    {
        if (options.array[i-1].processed)
            continue;

        optType ot = options.array[i-1].option;
        er val = options.array[i-1].value.get();
        if (ot == opt_style_string)
        {
            style->push_pair(ot, ecopy(val));
            options.array[i-1].processed = true;
        }
    }

    cellgrouping = style->get_cellgrouping(opt_CellGrouping);
    celllabelautodelete = style->get_bool(opt_CellLabelAutoDelete);
    cellautooverwrite = style->get_bool(opt_CellAutoOverwrite);
    celllabelmargin = style->get_rectangle(opt_CellLabelMargin);
    cellmargin = style->get_rectangle(opt_CellMargin);
    cellframe = style->get_rectangle(opt_CellFrame);
    cellframemargin = style->get_rectangle(opt_CellFrameMargin);
    cellframecolor = style->get_color(opt_CellFrameColor);

    for (ulong i = options.array.size(); i > 0 ; i--)
    {
        if (options.array[i-1].processed)
            continue;

        optType ot = options.array[i-1].option;
        er val = options.array[i-1].value.get();
        switch (ot)
        {
            case opt_CellGrouping:
                style->update_cellgrouping(cellgrouping, ot, val);
                break;
            case opt_CellLabelAutoDelete:
                style->update_bool(celllabelautodelete, ot, val);
                break;
            case opt_CellLabelMargin:
                style->update_rectangle(celllabelmargin, ot, val);
                break;
            case opt_CellAutoOverwrite:
                style->update_bool(cellautooverwrite, ot, val);
                break;
            case opt_CellMargin:
                style->update_rectangle(cellmargin, ot, val);
                break;
            case opt_CellFrame:
                style->update_rectangle(cellframe, ot, val);
                break;
            case opt_CellFrameMargin:
                style->update_rectangle(cellframemargin, ot, val);
                break;
            case opt_CellFrameColor:
                style->update_color(cellframecolor, ot, val);
                break;
            default:
                style->push_pair(ot, ecopy(val));
        }
    }

    style->pop_options();
    
    options.restore();

    flags |= BNFLAG_MEASURED;

//std::cout << "cellbox::measure_quick returning" << std::endl;
}


void cellbox::measure(boxmeasurearg ma)
{
    if (mexpr != nullptr)
    {
        //fcolor = da.style.parent_nb()->cCellForegroundInput;
        mexpr->measure(ma);
        bracket_offy = 30;
        bracket_sizey = mexpr->sizey;
        sizex = mexpr->sizex + 2*bracket_offy;
        sizey = mexpr->sizey + 2*bracket_offy;
        centery = 0.5*sizey;
        return;
    }


//std::cout << "cellbox::measure called" << std::endl;

    ma.level++;

//std::cout << "cellbox::measure options: " << ex_tostring_full(options.get()) << std::endl;
//std::cout << "cellbox::measure   style: " << ex_tostring_full(style.get()) << std::endl;

    notebook * nb = ma.style->parent_nb();

    if (label.cbox != nullptr)
    {
        delete label.cbox;
        label.cbox = nullptr;
    }

    options.set_restore_point();
    options.append_ex_list(ma.style->get_opt(opt_CellOptions));

    // measure label first

    wex cell_label_options(ecopy(ma.style->get_opt(opt_CellLabelOptions)));

    for (ulong i = options.array.size(); i > 0 ; i--)
    {
        optType ot = options.array[i-1].option;
        er val = options.array[i-1].value.get();
        if (ot == opt_CellLabel)
        {
            if (label.cbox != nullptr)
                delete label.cbox;
            label.cbox = rowbox_from_ex(val);
            options.array[i-1].processed = true;
        }
        else if (ot == opt_CellLabelOptions)
        {
            // TODO: check for "expired" cell
            if (ehas_head_sym_length(val, gs.sym_sList.get(), 2))
                val = echild(val,1);
            cell_label_options.reset(ecopy(val));
            options.array[i-1].processed = true;
        }
        else
        {
            options.array[i-1].processed = false;
        }
    }

    if (label.cbox != nullptr)
    {
        ma.style->start_push_options();
        ma.style->push_ex_list(cell_label_options.get());
        ma.style->finish_push_options();
        label.cbox->measure(ma);
        ma.style->pop_options();
    }

    // measure body

    ma.style->start_push_options();

    for (ulong i = options.array.size(); i > 0 ; i--)
    {
        if (options.array[i-1].processed)
            continue;

        optType ot = options.array[i-1].option;
        er val = options.array[i-1].value.get();
        if (ot == opt_style_string)
        {
            ma.style->push_pair(ot, ecopy(val));
            options.array[i-1].processed = true;
        }
    }

    cellgrouping = ma.style->get_cellgrouping(opt_CellGrouping);
    celllabelautodelete = ma.style->get_bool(opt_CellLabelAutoDelete);
    cellautooverwrite = ma.style->get_bool(opt_CellAutoOverwrite);
    celllabelmargin = ma.style->get_rectangle(opt_CellLabelMargin);
    cellmargin = ma.style->get_rectangle(opt_CellMargin);
    cellframe = ma.style->get_rectangle(opt_CellFrame);
    cellframemargin = ma.style->get_rectangle(opt_CellFrameMargin);
    cellframecolor = ma.style->get_color(opt_CellFrameColor);

    for (ulong i = options.array.size(); i > 0 ; i--)
    {
        if (options.array[i-1].processed)
            continue;

        optType ot = options.array[i-1].option;
        er val = options.array[i-1].value.get();
        switch (ot)
        {
            case opt_CellGrouping:
                ma.style->update_cellgrouping(cellgrouping, ot, val);
                break;
            case opt_CellLabelAutoDelete:
                ma.style->update_bool(celllabelautodelete, ot, val);
                break;
            case opt_CellLabelMargin:
                ma.style->update_rectangle(celllabelmargin, ot, val);
                break;
            case opt_CellAutoOverwrite:
                ma.style->update_bool(cellautooverwrite, ot, val);
                break;
            case opt_CellMargin:
                ma.style->update_rectangle(cellmargin, ot, val);
                break;
            case opt_CellFrame:
                ma.style->update_rectangle(cellframe, ot, val);
                break;
            case opt_CellFrameMargin:
                ma.style->update_rectangle(cellframemargin, ot, val);
                break;
            case opt_CellFrameColor:
                ma.style->update_color(cellframecolor, ot, val);
                break;
            default:
                ma.style->push_pair(ot, ecopy(val));
        }
    }

    ma.style->finish_push_options();

    cellmargin.scale(nb->magnification);
    cellframe.scale(nb->magnification);
    cellframemargin.scale(nb->magnification);
    celllabelmargin.scale(nb->magnification);

//std::cout << "cellmargin: "; cellmargin.print(); std::cout << std::endl;
//std::cout << "cellframe: "; cellframe.print(); std::cout << std::endl;
//std::cout << "cellframemargin: "; cellframemargin.print(); std::cout << std::endl;
//std::cout << "celllabelmargin: "; celllabelmargin.print(); std::cout << std::endl;

    double fs = ma.style->get_font();
    double cellbracket_w = glb_dingbat.get_char_sizex(DINGBAT_CELLGEN, fontsize_size(fs));
    double cellbracket_h = glb_dingbat.get_char_sizey(DINGBAT_CELLGEN, fontsize_size(fs));
//std::cout << "cellbox::measure cellbracket_w: " << cellbracket_w << std::endl;
//std::cout << "cellbox::measure cellbracket_h: " << cellbracket_h << std::endl;

    double pad_left = cellmargin.left + cellframe.left + cellframemargin.left;
    double pad_right = cellmargin.right + cellframe.right + cellframemargin.right;

    body.cbox->measure(boxmeasurearg(ma.style,
                                     ma.wrap_width - pad_left - pad_right - cellbracket_w*(ma.level + 3),
                                     ma.level));

    ma.style->pop_options();

    if (label.cbox != nullptr)
    {
        if (label.cbox->sizex + celllabelmargin.left + celllabelmargin.right <= cellmargin.left)
        {
            // label on left, align centers and extend cellmargin.top/bottom if needed
            label.offx = cellmargin.left - (label.cbox->sizex + celllabelmargin.right);

            // need celllabelmargin.top  <= label.offy
            //  and label.offy + celllabel.cbox->center =
            //      cellmargin.top + cellframe.top + cellframemargin.top + body.cbox->centery

            label.offy = cellmargin.top + cellframe.top + cellframemargin.top + body.cbox->centery - label.cbox->centery;
            if (label.offy < celllabelmargin.top)
            {
                cellmargin.top += celllabelmargin.top - label.offy;
                label.offy = celllabelmargin.top;
            }

            // need label.offy + label.cbox->sizey + celllabelmargin.bottom <=
            // cellmargin.top + cellframe.top + cellframemargin.top + body.cbox->sizey + cellmargin.bottom + cellframe.bottom + cellframemargin.bottom
            cellmargin.bottom = std::max(cellmargin.bottom, (label.offy + label.cbox->sizey + celllabelmargin.bottom) -
                                                            (cellmargin.top + cellframe.top + cellframemargin.top) -
                                                            body.cbox->sizey -
                                                            (cellframe.bottom + cellframemargin.bottom));
        }
        else
        {
            // label on top
            cellmargin.top = std::max(cellmargin.top, label.cbox->sizey + celllabelmargin.top + celllabelmargin.bottom);
            label.offy = cellmargin.top - (label.cbox->sizey + celllabelmargin.bottom);
            label.offx = celllabelmargin.left;
        }
    }

    double pad_top = cellmargin.top + cellframe.top + cellframemargin.top;
    double pad_bottom = cellmargin.bottom + cellframe.bottom + cellframemargin.bottom;

    body.offx = pad_left;
    body.offy = pad_top;
    sizex = pad_left + body.cbox->sizex + pad_right;
    sizey = pad_top + body.cbox->sizey + pad_bottom;
    centery = pad_top + body.cbox->centery;

    bracket_offy = cellmargin.top;
    bracket_sizey = sizey - cellmargin.bottom - cellmargin.top;

    options.restore();

    flags |= BNFLAG_MEASURED;

//std::cout << "cellbox::measure returning" << std::endl;
}



void cellbox::draw_pre(boxdrawarg da)
{
    if (mexpr != nullptr)
    {
        mexpr->draw_pre(boxdrawarg(da, bracket_offy, bracket_offy));
        return;
    }

    // TODO: draw frame and background

    body.cbox->draw_pre(boxdrawarg(da, body.offx, body.offy));        
}

void cellbox::draw_main(boxdrawarg da)
{
    _draw_cellgroup_bracket(this, da, false);

    if (mexpr != nullptr)
    {
        mexpr->draw_main(boxdrawarg(da, bracket_offy, bracket_offy));
        return;
    }

    da.dflags = (get_cellgrouping().first == cellgt_Input) ? DFLAG_SCOLOR : 0;
    body.cbox->draw_main(boxdrawarg(da, body.offx, body.offy, 0));
    da.dflags = 0;
    if (label.cbox != nullptr)
        label.cbox->draw_main(boxdrawarg(da, label.offx, label.offy, DFLAG_IGNORESEL));
}

void cellbox::draw_post(boxdrawarg da)
{
    if (mexpr != nullptr)
    {
        mexpr->draw_post(boxdrawarg(da, 0, 0));
        return;
    }

    da.dflags = (get_cellgrouping().first == cellgt_Input) ? DFLAG_SCOLOR : 0;
    body.cbox->draw_post(boxdrawarg(da, body.offx, body.offy));        
}
