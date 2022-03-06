#include "box_convert.h"
#include "notebook.h"
#include "ex_parse.h"
#include "ex_cont.h"
#include "eval.h"
#include "arithmetic.h"
#include "ex_print.h"
#include "timing.h"

static void gui_error(const std::string & s)
{
//    std::cerr << "<!gui>: " << s << std::endl; 
//    SleepMS(100);
}


void parsedoptions::remove(optType o)
{
    for (ulong i = 0; i < array.size(); i++)
    {
        if (array[i].option == o)
        {
            std::swap(array[i], array.back());
            array.pop_back();
            return;
        }
    }
};

void parsedoptions::push_ex_rule(er e)
{
    if (ehas_head_sym_length(e, gs.sym_sRule.get(), 2))
    {
        er a = echild(e,1);
        if (eis_str(a, "FontSize") || eis_sym(a, gs.sym_sFontSize.get()))
        {
            push_pair(opt_FontSize, ecopychild(e,2));
            return;
        }
        else if (eis_str(a, "FontFamily") || eis_sym(a, gs.sym_sFontFamily.get()))
        {
            push_pair(opt_FontFamily, ecopychild(e,2));
            return;
        }
        else if (eis_str(a, "FontWeight") || eis_sym(a, gs.sym_sFontWeight.get()))
        {
            push_pair(opt_FontWeight, ecopychild(e,2));
            return;
        }
        else if (eis_str(a, "FontSlant") || eis_sym(a, gs.sym_sFontSlant.get()))
        {
            push_pair(opt_FontSlant, ecopychild(e,2));
            return;
        }
        else if (eis_str(a, "FontColor") || eis_sym(a, gs.sym_sFontColor.get()))
        {
            push_pair(opt_FontColor, ecopychild(e,2));
            return;
        }
        else if (eis_str(a, "FontOpacity")/* || eis_sym(a, gs.sym_sFontOpacity.get())*/)
        {
            push_pair(opt_FontOpacity, ecopychild(e,2));
            return;
        }
        else if (eis_str(a, "AutoSpacing") || eis_sym(a, gs.sym_sAutoSpacing.get()))
        {
            push_pair(opt_AutoSpacing, ecopychild(e,2));
            return;
        }
        else if (eis_str(a, "CellLabelMargin")/* || eis_sym(a, gs.sym_sCellLabelMargin.get())*/)
        {
            push_pair(opt_CellLabelMargin, ecopychild(e,2));
            return;
        }
        else if (eis_str(a, "CellLabelAutoDelete")/* || eis_sym(a, gs.sym_sCellLabelAutoDelete.get())*/)
        {
            push_pair(opt_CellLabelAutoDelete, ecopychild(e,2));
            return;
        }
        else if (eis_str(a, "CellAutoOverwrite")/* || eis_sym(a, gs.sym_sCellAutoOverwrite.get())*/)
        {
            push_pair(opt_CellAutoOverwrite, ecopychild(e,2));
            return;
        }
        else if (eis_str(a, "CellGrouping")/* || eis_sym(a, gs.sym_sCellGrouping.get())*/)
        {
            push_pair(opt_CellGrouping, ecopychild(e,2));
            return;
        }
        else if (eis_str(a, "CellMargin")/* || eis_sym(a, gs.sym_sCellMargin.get())*/)
        {
            push_pair(opt_CellMargin, ecopychild(e,2));
            return;
        }
        else if (eis_str(a, "CellFrame")/* || eis_sym(a, gs.sym_sCellFrame.get())*/)
        {
            push_pair(opt_CellFrame, ecopychild(e,2));
            return;
        }
        else if (eis_str(a, "CellFrameMargin")/* || eis_sym(a, gs.sym_sCellFrameMargin.get())*/)
        {
            push_pair(opt_CellFrameMargin, ecopychild(e,2));
            return;
        }
        else if (eis_str(a, "CellFrameColor")/* || eis_sym(a, gs.sym_sCellFrameColor.get())*/)
        {
            push_pair(opt_CellFrameColor, ecopychild(e,2));
            return;
        }
        else if (eis_str(a, "CellBraketColor")/* || eis_sym(a, gs.sym_sCellBraketColor.get())*/)
        {
            push_pair(opt_CellBracketColor, ecopychild(e,2));
            return;
        }
        else if (eis_str(a, "CellLabel") || eis_sym(a, gs.sym_sCellLabel.get()))
        {
            push_pair(opt_CellLabel, ecopychild(e,2));
            return;
        }
        else if (eis_str(a, "CellOptions")/* || eis_sym(a, gs.sym_sCellOptions.get())*/)
        {
            push_pair(opt_CellOptions, ecopychild(e,2));
            return;
        }
        else if (eis_str(a, "CellLabelOptions")/* || eis_sym(a, gs.sym_sCellLabelOptions.get())*/)
        {
            push_pair(opt_CellLabelOptions, ecopychild(e,2));
            return;
        }
        else if (eis_str(a, "GridBoxOptions")/* || eis_sym(a, gs.sym_sGridBoxOptions.get())*/)
        {
            push_pair(opt_GridBoxOptions, ecopychild(e,2));
            return;
        }
    }
    else if (ehas_head_sym(e, gs.sym_sList.get()))
    {
        for (ulong i = elength(e); i > 0; i--)
            push_ex_rule(echild(e,i-1));
        return;
    }
    else if (eis_str(e))
    {
        push_pair(opt_style_string, ecopy(e));
        return;
    }

    gui_error("parsedoptions::push_ex_rule: improper option " + ex_tostring_full(e));
}

void parsedoptions::append_ex_list(er e)
{
    // TODO: also can use Directive
    if (!ehas_head_sym(e, gs.sym_sList.get()))
        return;

    for (ulong i = elength(e); i > 0; i--)
        push_ex_rule(echild(e,i-1));
}

void parsedoptions::get_ex(std::vector<wex> & l)
{
//std::cout << "parsedoptions::get_ex called" << std::endl;

    for (ulong i = array.size(); i > 0; i--)
    {
        er v = array[i-1].value.get();


//std::cout << "option: " << array[i-1].option << std::endl;

        switch (array[i-1].option)
        {
            case opt_FontSize:
                l.push_back(wex(emake_node(gs.sym_sRule.copy(), emake_str("FontSize"), ecopy(v))));
                break;
            case opt_FontFamily:
                l.push_back(wex(emake_node(gs.sym_sRule.copy(), emake_str("FontFamily"), ecopy(v))));
                break;
            case opt_FontWeight:
                l.push_back(wex(emake_node(gs.sym_sRule.copy(), emake_str("FontWeight"), ecopy(v))));
                break;
            case opt_FontSlant:
                l.push_back(wex(emake_node(gs.sym_sRule.copy(), emake_str("FontSlant"), ecopy(v))));
                break;
            case opt_FontColor:
                l.push_back(wex(emake_node(gs.sym_sRule.copy(), emake_str("FontColor"), ecopy(v))));
                break;
            case opt_FontOpacity:
                l.push_back(wex(emake_node(gs.sym_sRule.copy(), emake_str("FontOpacity"), ecopy(v))));
                break;
            case opt_AutoSpacing:
                l.push_back(wex(emake_node(gs.sym_sRule.copy(), emake_str("AutoSpacing"), ecopy(v))));
                break;
            case opt_CellLabelMargin:
                l.push_back(wex(emake_node(gs.sym_sRule.copy(), emake_str("CellLabelMargin"), ecopy(v))));
                break;
            case opt_CellLabelAutoDelete:
                l.push_back(wex(emake_node(gs.sym_sRule.copy(), emake_str("CellLabelAutoDelete"), ecopy(v))));
                break;
            case opt_CellAutoOverwrite:
                l.push_back(wex(emake_node(gs.sym_sRule.copy(), emake_str("CellAutoOverwrite"), ecopy(v))));
                break;
            case opt_CellGrouping:
                l.push_back(wex(emake_node(gs.sym_sRule.copy(), emake_str("CellGrouping"), ecopy(v))));
                break;
            case opt_CellMargin:
                l.push_back(wex(emake_node(gs.sym_sRule.copy(), emake_str("CellMargin"), ecopy(v))));
                break;
            case opt_CellFrame:
                l.push_back(wex(emake_node(gs.sym_sRule.copy(), emake_str("CellFrame"), ecopy(v))));
                break;
            case opt_CellFrameMargin:
                l.push_back(wex(emake_node(gs.sym_sRule.copy(), emake_str("CellFrameMargin"), ecopy(v))));
                break;
            case opt_CellFrameColor:
                l.push_back(wex(emake_node(gs.sym_sRule.copy(), emake_str("CellFrameColor"), ecopy(v))));
                break;
            case opt_CellBracketColor:
                l.push_back(wex(emake_node(gs.sym_sRule.copy(), emake_str("CellBracketColor"), ecopy(v))));
                break;
            case opt_CellLabel:
                l.push_back(wex(emake_node(gs.sym_sRule.copy(), emake_str("CellLabel"), ecopy(v))));
                break;
            case opt_style_string:
                l.push_back(wex(ecopy(v)));
                break;
            case opt_CellOptions:
                l.push_back(wex(emake_node(gs.sym_sRule.copy(), emake_str("CellOptions"), ecopy(v))));
                break;
            case opt_CellLabelOptions:
                l.push_back(wex(emake_node(gs.sym_sRule.copy(), emake_str("CellLabelOptions"), ecopy(v))));
                break;
            case opt_GridBoxOptions:
                l.push_back(wex(emake_node(gs.sym_sRule.copy(), emake_str("GridBoxOptions"), ecopy(v))));
                break;
            default:
                std::cerr << "unhandled case in parsedoptions.get_ex" << std::endl;
                assert(false);
        }
    }
}




stylestack::stylestack(notebook * nb_) : nb(nb_), const1_ok(false)
{
    for (int i = 0; i < opt_MAX; i++)
        opts_idxs[i] = -1;
}


ex etry_double(er e)
{
    if (eis_double(e))
        return ecopy(e);

    double d = econvert_todouble(e);
    return std::isfinite(d) ? emake_double(d) : nullptr;
}

ex etry_color(er e)
{
    if (!ehas_head_sym_length(e, gs.sym_sRGBColor.get(), 3))
        return nullptr;

    if (eis_double(echild(e,1)) && eis_double(echild(e,2)) && eis_double(echild(e,3)))
        return ecopy(e);

    uex r(etry_double(echild(e,1)));
    uex g(etry_double(echild(e,2)));
    uex b(etry_double(echild(e,3)));

    if (r.get() == nullptr || g.get() == nullptr || b.get() == nullptr)
        return nullptr;

    return emake_node(gs.sym_sRGBColor.copy(), r.release(), g.release(), b.release());
}

bool stylestack::get_bool(optType opt)
{
//std::cout << "stylestack::get_bool opt = " << opt << std::endl;
//std::cout << "                     ans = " << eis_sym(get_opt(opt), gs.sym_sTrue.get()) << std::endl;


    return eis_sym(get_opt(opt), gs.sym_sTrue.get());
}

uint32_t stylestack::get_color(optType opt)
{
    er e = get_opt(opt);

    if (!ehas_head_sym_length(e, gs.sym_sRGBColor.get(), 3))
        return 0x80808080;

    uex R(etry_double(echild(e,1)));
    uex G(etry_double(echild(e,2)));
    uex B(etry_double(echild(e,3)));

    if (R.get() == nullptr || G.get() == nullptr || B.get() == nullptr)
        return 0x808080;

    uint32_t r = 255.99*edouble_number(R.get());
    uint32_t g = 255.99*edouble_number(G.get());
    uint32_t b = 255.99*edouble_number(B.get());
    r = std::min(r, uint32_t(255));
    g = std::min(g, uint32_t(255));
    b = std::min(b, uint32_t(255));
    return RGB_COLOR(r, g, b);
}



ex eunpack(er e, ulong depth, ulong maxdepth)
{
	if (depth > maxdepth)
		return ecopy(e);

	uex r(ecopy(e));
	if (eis_parray(e))
		r.setnz(ecopy(eparray_get_normal(e)));

	if (eis_node(r.get()))
	{
		ulong n = elength(r.get());
		uex f; f.init_push_backx(eunpack(r.child(0), depth + 1, maxdepth), n);
		for (ulong i = 1; i <= n; i++)
			f.push_back(eunpack(r.child(i), depth + 1, maxdepth));
		r.setnz(f.release());
	}

	return r.release();
}


double econvert_to_double(er x, double fail)
{
    double ret = fail;
    switch (etype(x))
    {
        case ETYPE_DOUBLE:
            ret = eto_double(x)->number;
            break;
        case ETYPE_INT:
            ret = fmpz_get_d(eint_data(x));
            break;
        case ETYPE_RAT:
            ret = fmpq_get_d(erat_data(x));
            break;
        case ETYPE_REAL:
            ret = arf_get_d(arb_midref(ereal_data(x)), ARF_RND_NEAR);
            break;
    }
    return std::isfinite(ret) ? ret : fail;
}


bool eis_rectangle(er e)
{
    if (ehas_head_sym_length(e, gs.sym_sList.get(), 2))
    {
        if (ehas_head_sym_length(echild(e,1), gs.sym_sList.get(), 2))
        {
            if (!std::isfinite(econvert_to_double(echild(e,1,1), NAN)))
                return false;
            if (!std::isfinite(econvert_to_double(echild(e,1,2), NAN)))
                return false;
        }
        else
        {
            return false;
        }

        if (ehas_head_sym_length(echild(e,2), gs.sym_sList.get(), 2))
        {
            if (!std::isfinite(econvert_to_double(echild(e,2,1), NAN)))
                return false;
            if (!std::isfinite(econvert_to_double(echild(e,2,2), NAN)))
                return false;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }

    return true;
}

rectangle stylestack::get_rectangle(optType opt)
{
    uex E(ecopy(get_opt(opt)));

    if (eis_parray(E.get()))
        E.setnz(eunpack(E.get(), 1, 2));

    er e = E.get();

//std::cout << "stylestack::get_rectangle " << ex_tostring_full(e) << std::endl;


    double left = 0.0;
    double right = 0.0;
    double bottom = 0.0;
    double top = 0.0;

    if (ehas_head_sym_length(e, gs.sym_sList.get(), 2))
    {
        if (ehas_head_sym_length(echild(e,1), gs.sym_sList.get(), 2))
        {
            left = econvert_to_double(echild(e,1,1), left);
            right = econvert_to_double(echild(e,1,2), right);
        }
        if (ehas_head_sym_length(echild(e,2), gs.sym_sList.get(), 2))
        {
            bottom = econvert_to_double(echild(e,2,1), bottom);
            top = econvert_to_double(echild(e,2,2), top);
        }
    }

//std::cout << "stylestack::get_rectangle " << left << ", " << right << ", " << bottom << ", " << top << std::endl;

    return {left, right, bottom, top};
}

std::pair<cellGroupType, int32_t> stylestack::get_cellgrouping(optType opt)
{
    er e = get_opt(opt);

    if (eis_str(e, "Input"))
    {
        return {cellgt_Input, 0};
    }
    else if (eis_str(e, "Output"))
    {
        return {cellgt_Output, 0};
    }
    else if (eis_str(e, "Message"))
    {
        return {cellgt_Message, 0};
    }
    else if (eis_str(e, "Print"))
    {
        return {cellgt_Print, 0};
    }
    else if (eis_str(e, "Text"))
    {
        return {cellgt_Text, 0};
    }
    else if (eis_str(e, "Title"))
    {
        return {cellgt_Title, 0};
    }
    else if (eis_str(e, "Subsubsection"))
    {
        return {cellgt_Subsubsection, 0};
    }
    else if (eis_str(e, "Subsection"))
    {
        return {cellgt_Subsection, 0};
    }
    else if (eis_str(e, "Section"))
    {
        return {cellgt_Section, 0};
    }
    else
    {
        return {cellgt_Default, 0};
    }
}


ex emake_color(double r, double g, double b);

void stylestack::push_pair(optType opt, ex V)
{
    uex VV(V);
    er v = VV.get();

    switch (opt)
    {
        case opt_FontSize:
        {
            ex t = etry_double(v);
            if (t != nullptr)
                push_opt(opt_FontSize, t);
            else
                gui_error("improper FontSize " + ex_tostring_full(v));
            break;
        }
        case opt_FontFamily:
        {
            if (eis_str(v))
                push_opt(opt_FontFamily, ecopy(v));
            else
                gui_error("improper FontFamily " + ex_tostring_full(v));
            break;
        }
        case opt_FontWeight:
        {
            if (eis_str(v))
                push_opt(opt_FontWeight, ecopy(v));
            else
                gui_error("improper FontWeight " + ex_tostring_full(v));
            break;
        }
        case opt_FontSlant:
        {
            if (eis_str(v))
                push_opt(opt_FontSlant, ecopy(v));
            else
                gui_error("improper FontSlant " + ex_tostring_full(v));
            break;
        }
        case opt_FontColor:
        {
            ex t = etry_color(v);
            if (t != nullptr)
                push_opt(opt_FontColor, t);
            else
                gui_error("improper FontColor " + ex_tostring_full(v));
            break;
        }
        case opt_AutoSpacing:
        {
            if (eis_sym(v, gs.sym_sTrue.get()) || eis_sym(v, gs.sym_sFalse.get()))
                push_opt(opt_AutoSpacing, ecopy(v));
            else
                gui_error("improper AutoSpacing " + ex_tostring_full(v));
            break;
        }
        case opt_CellLabelAutoDelete:
        {
            if (eis_sym(v, gs.sym_sTrue.get()) || eis_sym(v, gs.sym_sFalse.get()))
                push_opt(opt_CellLabelAutoDelete, ecopy(v));
            else
                gui_error("improper CellLabelAutoDelete " + ex_tostring_full(v));
            break;
        }
        case opt_CellAutoOverwrite:
        {
            if (eis_sym(v, gs.sym_sTrue.get()) || eis_sym(v, gs.sym_sFalse.get()))
                push_opt(opt_CellAutoOverwrite, ecopy(v));
            else
                gui_error("improper CellAutoOverwrite " + ex_tostring_full(v));
            break;
        }
        case opt_CellGrouping:
        {
            if (eis_str(v))
                push_opt(opt_CellGrouping, ecopy(v));
            else
                gui_error("improper CellGrouping " + ex_tostring_full(v));
            break;
        }
        case opt_CellLabelMargin:
        case opt_CellMargin:
        case opt_CellFrame:
        case opt_CellFrameMargin:
        {

            if (eis_rectangle(v))
                push_opt(opt, ecopy(v));
            else
                gui_error("improper rectangle " + ex_tostring_full(v));
            break;
        }
        case opt_style_string:
        {
            // TODO parse these from somewhere else
            if (eis_str(v, "Input"))
            {
                push_pair(opt_FontSize, emake_double(11.0));
                push_pair(opt_FontFamily, emake_str("Courier"));
                push_pair(opt_FontWeight, emake_str("Bold"));
                push_pair(opt_FontSlant, emake_str("Plain"));
                push_pair(opt_FontColor, gs.const_color_black.copy());
                push_pair(opt_AutoSpacing, gs.sym_sTrue.copy());
                push_pair(opt_CellLabelAutoDelete, gs.sym_sTrue.copy());
                push_pair(opt_CellAutoOverwrite, gs.sym_sFalse.copy());
                push_pair(opt_CellGrouping, emake_str("Input"));
                push_pair(opt_CellMargin, emake_list(emake_list(emake_double(72.0), emake_double(5.0)),
                                                     emake_list(emake_double(10.0), emake_double(8.0))));
            }
            else if (eis_str(v, "Output"))
            {
                push_pair(opt_FontSize, emake_double(11.0));
                push_pair(opt_FontFamily, emake_str("Courier"));
                push_pair(opt_FontWeight, emake_str("Regular"));
                push_pair(opt_FontColor, gs.const_color_black.copy());
                push_pair(opt_AutoSpacing, gs.sym_sTrue.copy());
                push_pair(opt_CellLabelAutoDelete, gs.sym_sTrue.copy());
                push_pair(opt_CellAutoOverwrite, gs.sym_sTrue.copy());
                push_pair(opt_CellGrouping, emake_str("Output"));
                push_pair(opt_CellMargin, emake_list(emake_list(emake_double(72.0), emake_double(5.0)),
                                                     emake_list(emake_double(12.0), emake_double(8.0))));
            }
            else if (eis_str(v, "Print"))
            {
                push_pair(opt_FontSize, emake_double(10.0));
                push_pair(opt_FontFamily, emake_str("Courier"));
                push_pair(opt_FontWeight, emake_str("Regular"));
                push_pair(opt_FontColor, gs.const_color_black.copy());
                push_pair(opt_AutoSpacing, gs.sym_sTrue.copy());
                push_pair(opt_CellAutoOverwrite, gs.sym_sTrue.copy());
                push_pair(opt_CellGrouping, emake_str("Print"));
                push_pair(opt_CellMargin, emake_list(emake_list(emake_double(72.0), emake_double(5.0)),
                                                     emake_list(emake_double(4.0), emake_double(4.0))));
            }
            else if (eis_str(v, "Title"))
            {
                push_pair(opt_FontSize, emake_double(32.0));
                push_pair(opt_FontFamily, emake_str("Segoe"));
                push_pair(opt_FontColor, emake_color(0xcc/255.0,0x0a/255.0,0x02/255.0));
                push_pair(opt_CellGrouping, emake_str("Title"));
                push_pair(opt_CellMargin, emake_list(emake_list(emake_double(16.0), emake_double(16.0)),
                                                     emake_list(emake_double(16.0), emake_double(16.0))));
            }
            else if (eis_str(v, "Section"))
            {
                push_pair(opt_FontSize, emake_double(24.0));
                push_pair(opt_FontFamily, emake_str("Segoe"));
                push_pair(opt_FontColor, emake_color(0xC2/255.0,0x4B/255.0,0x15/255.0));
                push_pair(opt_CellGrouping, emake_str("Section"));
                push_pair(opt_CellMargin, emake_list(emake_list(emake_double(32.0), emake_double(16.0)),
                                                     emake_list(emake_double(14.0), emake_double(14.0))));
            }
            else if (eis_str(v, "Subsection"))
            {
                push_pair(opt_FontSize, emake_double(20.0));
                push_pair(opt_FontFamily, emake_str("Segoe"));
                push_pair(opt_FontColor, emake_color(0xC6/255.0,0x6B/255.0,0x29/255.0));
                push_pair(opt_CellGrouping, emake_str("Subsection"));
                push_pair(opt_CellMargin, emake_list(emake_list(emake_double(48.0), emake_double(16.0)),
                                                     emake_list(emake_double(12.0), emake_double(12.0))));
            }
            else if (eis_str(v, "Subsubsection"))
            {
                push_pair(opt_FontSize, emake_double(16.0));
                push_pair(opt_FontFamily, emake_str("Segoe"));
                push_pair(opt_FontColor, emake_color(0xB6/255.0,0x37/255.0,0x08/255.0));
                push_pair(opt_CellGrouping, emake_str("Subsubsection"));
                push_pair(opt_CellMargin, emake_list(emake_list(emake_double(64.0), emake_double(16.0)),
                                                     emake_list(emake_double(10.0), emake_double(10.0))));
            }
            else if (eis_str(v, "Text"))
            {
                push_pair(opt_FontSize, emake_double(11.0));
                push_pair(opt_FontFamily, emake_str("Verdana"));
                push_pair(opt_CellGrouping, emake_str("Text"));
                push_pair(opt_CellMargin, emake_list(emake_list(emake_double(32.0), emake_double(16.0)),
                                                     emake_list(emake_double(8.0), emake_double(8.0))));
            }
            else if (eis_str(v, "ExampleText"))
            {
                push_pair(opt_FontSize, emake_double(11.0));
                push_pair(opt_FontFamily, emake_str("Verdana"));
                push_pair(opt_CellGrouping, emake_str("Text"));
                push_pair(opt_CellMargin, emake_list(emake_list(emake_double(72.0), emake_double(16.0)),
                                                     emake_list(emake_double(8.0), emake_double(8.0))));
            }
            else if (eis_str(v, "Message"))
            {
                push_pair(opt_FontSize, emake_double(11.0));
                push_pair(opt_FontFamily, emake_str("Segoe"));
                push_pair(opt_FontColor, emake_color(0xb0/255.0,0x50/255.0,0x18/255.0));
                push_pair(opt_CellAutoOverwrite, gs.sym_sTrue.copy());
                push_pair(opt_CellGrouping, emake_str("Message"));
                push_pair(opt_CellMargin, emake_list(emake_list(emake_double(72.0), emake_double(16.0)),
                                                     emake_list(emake_double(8.0), emake_double(8.0))));
            }
            else if (eis_str(v, "CellLabel"))
            {
                push_pair(opt_FontSize, emake_double(8.0));
                push_pair(opt_FontFamily, emake_str("Segoe"));
                push_pair(opt_FontColor, emake_color(0.25,7.0/16,0.50));
            }
            else
            {
                gui_error("unhandled style " + ex_tostring_full(v));
            }
            break;
        }
        default:
        {
            gui_error("unhandled option " + stdstring_tostring(opt) + " -> " + ex_tostring_full(v));
            break;
        }
    }
}

void stylestack::push_ex_rule(er e)
{
    if (ehas_head_sym_length(e, gs.sym_sRule.get(), 2))
    {
        er a = echild(e,1);
        er b = echild(e,2);
        if (eis_str(a, "FontSize") || eis_sym(a, gs.sym_sFontSize.get()))
        {
            push_pair(opt_FontSize, ecopy(b));
            return;
        }
        else if (eis_str(a, "FontFamily") || eis_sym(a, gs.sym_sFontFamily.get()))
        {
            push_pair(opt_FontFamily, ecopy(b));
            return;
        }
        else if (eis_str(a, "FontWeight") || eis_sym(a, gs.sym_sFontWeight.get()))
        {
            push_pair(opt_FontWeight, ecopy(b));
            return;
        }
        else if (eis_str(a, "FontSlant") || eis_sym(a, gs.sym_sFontSlant.get()))
        {
            push_pair(opt_FontSlant, ecopy(b));
            return;
        }
        else if (eis_str(a, "FontColor") || eis_sym(a, gs.sym_sFontColor.get()))
        {
            push_pair(opt_FontColor, ecopy(b));
            return;
        }
        else if (eis_str(a, "FontOpacity")/* || eis_sym(a, gs.sym_sFontOpacity.get())*/)
        {
            push_pair(opt_FontOpacity, ecopy(b));
            return;
        }
        else if (eis_str(a, "AutoSpacing") || eis_sym(a, gs.sym_sAutoSpacing.get()))
        {
            push_pair(opt_AutoSpacing, ecopy(b));
            return;
        }
        else if (eis_str(a, "CellLabelMargin")/* || eis_sym(a, gs.sym_sCellLabelMargin.get())*/)
        {
            push_pair(opt_CellLabelMargin, ecopy(b));
            return;
        }
        else if (eis_str(a, "CellLabelAutoDelete")/* || eis_sym(a, gs.sym_sCellLabelAutoDelete.get())*/)
        {
            push_pair(opt_CellLabelAutoDelete, ecopy(b));
            return;
        }
        else if (eis_str(a, "CellAutoOverwrite")/* || eis_sym(a, gs.sym_sCellAutoOverwrite.get())*/)
        {
            push_pair(opt_CellAutoOverwrite, ecopy(b));
            return;
        }
        else if (eis_str(a, "CellGrouping")/* || eis_sym(a, gs.sym_sCellGrouping.get())*/)
        {
            push_pair(opt_CellGrouping, ecopy(b));
            return;
        }
        else if (eis_str(a, "CellMargin")/* || eis_sym(a, gs.sym_sCellMargin.get())*/)
        {
            push_pair(opt_CellMargin, ecopy(b));
            return;
        }
        else if (eis_str(a, "CellFrame")/* || eis_sym(a, gs.sym_sCellFrame.get())*/)
        {
            push_pair(opt_CellFrame, ecopy(b));
            return;
        }
        else if (eis_str(a, "CellFrameMargin")/* || eis_sym(a, gs.sym_sCellFrameMargin.get())*/)
        {
            push_pair(opt_CellFrameMargin, ecopy(b));
            return;
        }
        else if (eis_str(a, "CellFrameColor")/* || eis_sym(a, gs.sym_sCellFrameColor.get())*/)
        {
            push_pair(opt_CellFrameColor, ecopy(b));
            return;
        }
        else if (eis_str(a, "CellBracketColor")/* || eis_sym(a, gs.sym_sCellBracketColor.get())*/)
        {
            push_pair(opt_CellBracketColor, ecopy(b));
            return;
        }
        else if (eis_str(a, "CellLabel")/* || eis_sym(a, gs.sym_sCellLabel.get())*/)
        {
            push_pair(opt_CellLabel, ecopy(b));
            return;
        }
        else if (eis_str(a, "CellOptions")/* || eis_sym(a, gs.sym_sCellOptions.get())*/)
        {
            push_pair(opt_CellOptions, ecopy(b));
            return;
        }
        else if (eis_str(a, "CellLabelOptions")/* || eis_sym(a, gs.sym_sCellLabelOptions.get())*/)
        {
            push_pair(opt_CellLabelOptions, ecopy(b));
            return;
        }
        else if (eis_str(a, "GridBoxOptions")/* || eis_sym(a, gs.sym_sGridBoxOptions.get())*/)
        {
            push_pair(opt_GridBoxOptions, ecopy(b));
            return;
        }
    }
    else if (ehas_head_sym(e, gs.sym_sList.get()))
    {
        for (ulong i = elength(e); i > 0; i--)
            push_ex_rule(echild(e,i-1));
        return;
    }
    else if (eis_str(e))
    {
        push_pair(opt_style_string, ecopy(e));
        return;
    }

    gui_error("stylestack::push_ex_rule: improper option " + ex_tostring_full(e));
}

void stylestack::push_ex_list(er e)
{
    new_options();
    if (!ehas_head_sym(e, gs.sym_sList.get()))
        return;

    for (ulong i = 0; i < elength(e); i++)
        push_ex_rule(echild(e,i+1));
}

void stylestack::pop_options()
{
    assert(!stack_dividers.empty());
    size_t n = stack_dividers.back();
    stack_dividers.pop_back();
    size_t i = value_stack.size();
    if (i > n)
    {
        do {
            i--;
            optType o = std::get<0>(value_stack[i]);
            opts_idxs[o] = std::get<1>(value_stack[i]);
            value_stack.pop_back();
        } while (i > n);
        init_const1();
    }
}

void stylestack::init_const1()
{
//std::cout << "stylestack::init_const1 called" << std::endl;
//print();

    er efamily = get_opt(opt_FontFamily);

    uint32_t family = FONTFAMILY_COURIER_REG;

    if (eis_str(efamily, "Courier"))
        family = FONTFAMILY_COURIER_REG;
    else if (eis_str(efamily, "Lucida"))
        family = FONTFAMILY_LUCIDA_REG;
    else if (eis_str(efamily, "Segoe"))
        family = FONTFAMILY_SEGOE_REG;
    else if (eis_str(efamily, "Times"))
        family = FONTFAMILY_TIMES_REG;
    else if (eis_str(efamily, "Verdana"))
        family = FONTFAMILY_VERDANA_REG;

    er eweight = get_opt(opt_FontWeight);
    if (eis_str(eweight, "Bold"))
        family &= ~(uint32_t)(1);

    const1_fontsize = edouble_number(get_opt(opt_FontSize));
    int32_t size = 0.5 + const1_fontsize*nb->magnification;
    size = std::max(size, MIN_FONT_SIZE);
    size = std::min(size, MAX_FONT_SIZE);

    const1_font = fontsize_make(family, size);

    const1_default_sizex = fontsize_default_sizex(const1_font);
    const1_default_sizey = fontsize_default_sizey(const1_font);
    const1_default_centery = fontsize_default_centery(const1_font);

    er e = get_opt(opt_FontColor);
    assert(ehas_head_sym_length(e, gs.sym_sRGBColor.get(), 3));

    uint32_t r = 255.99*edouble_number(echild(e,1));
    uint32_t g = 255.99*edouble_number(echild(e,2));
    uint32_t b = 255.99*edouble_number(echild(e,3));
    r = std::min(r, uint32_t(255));
    g = std::min(g, uint32_t(255));
    b = std::min(b, uint32_t(255));
    const1_fontcolor = RGB_COLOR(r, g, b);


    const1_ok = true;

//printf("     const1_font: %08lx\n", const1_font);
//printf("const1_fontcolor: %08lx\n", const1_fontcolor);
//printf("    const1_sizes: %f, %f, %f\n", const1_default_sizex, const1_default_sizey, const1_default_centery);
}

bool stylestack::get_autospacing()
{
//std::cout << "get_autospacing called" << std::endl;
    bool r = eis_sym(get_opt(opt_AutoSpacing), gs.sym_sTrue.get());
//std::cout << "get_autospacing returning " << r << std::endl;
    return r;
}

void stylestack::print() const
{
    std::cout << "opts_idx: ";
    for (int i = 0; i < opt_MAX; i++)
        std::cout << opts_idxs[i] << ", ";
    std::cout << std::endl;

    for (size_t i = 0; i < value_stack.size(); i++)
    {
        std::cout << "[" << i << "]:";
        std::cout << " opt: " << int(std::get<0>(value_stack[i]));
        std::cout << " prev: " << std::get<1>(value_stack[i]);
        std::cout << " value: " << ex_tostring_full(std::get<2>(value_stack[i]).get()) << std::endl;
    }

    std::cout << "dividers: ";
    for (size_t i = 0; i < stack_dividers.size(); i++)
        std::cout << stack_dividers[i] << ", ";
    std::cout << std::endl;
}
