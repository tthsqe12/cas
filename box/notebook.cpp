#include "font.h"
#include "graphics.h"
#include "notebook.h"
#include "box_convert.h"
#include "ex_print.h"

ex emake_color(double r, double g, double b)
{
    return emake_node(gs.sym_sRGBColor.copy(), emake_double(r), emake_double(g), emake_double(b));
}

notebook::notebook()
    : default_style(this)
{
    root = new rootbox(this, 0, 0, 0);
    copy_paste_buffer = nullptr;
    cursor_needs_fitting = true;
    offx = 0;
    offy = 0;
    window_sizex = 100;
    window_sizey = 100;
    filestring.clear();

    set_zoom128(128);

    default_style.push_opt(opt_FontSize, emake_double(11.0));
    default_style.push_opt(opt_FontFamily, emake_str("Times"));
    default_style.push_opt(opt_FontWeight, emake_str("Regular"));
    default_style.push_opt(opt_FontSlant, emake_str("Plain"));
    default_style.push_opt(opt_FontColor, gs.const_color_black.copy());
    default_style.push_opt(opt_FontOpacity, gs.const_double_one.copy());
    default_style.push_opt(opt_AutoSpacing, gs.sym_sFalse.copy());
    default_style.push_opt(opt_CellLabelAutoDelete, gs.sym_sFalse.copy());
    default_style.push_opt(opt_CellAutoOverwrite, gs.sym_sFalse.copy());
    default_style.push_opt(opt_CellGrouping, emake_node(gs.sym_sList.copy(), emake_str("text"), emake_cint(0)));
    default_style.push_opt(opt_CellLabelMargin, emake_list(emake_list(emake_double(6.0), emake_double(6.0)),
                                                           emake_list(emake_double(2.0), emake_double(2.0))));
    default_style.push_opt(opt_CellMargin, emake_list(emake_list(gs.const_double_zero.copy(), gs.const_double_zero.copy()), emake_list(gs.const_double_zero.copy(), gs.const_double_zero.copy())));
    default_style.push_opt(opt_CellFrame,emake_list(emake_list(gs.const_double_zero.copy(), gs.const_double_zero.copy()), emake_list(gs.const_double_zero.copy(), gs.const_double_zero.copy())));
    default_style.push_opt(opt_CellFrameMargin, emake_list(emake_list(gs.const_double_zero.copy(), gs.const_double_zero.copy()), emake_list(gs.const_double_zero.copy(), gs.const_double_zero.copy())));
    default_style.push_opt(opt_CellFrameColor, gs.const_color_black.copy());
    default_style.push_opt(opt_CellBracketColor, emake_node(gs.sym_sRGBColor.copy(), emake_double(0x90/255.0), emake_double(0x80/255.0), emake_double(0xc0/255.0)));
    default_style.push_opt(opt_CellLabel, gs.sym_sNone.copy());

    default_style.push_opt(opt_style_string, emake_str("Input"));

    default_style.push_opt(opt_CellOptions, gs.const_empty_list.copy());
    default_style.push_opt(opt_CellLabelOptions, emake_list(emake_str("CellLabel")));
    default_style.push_opt(opt_GridBoxOptions, gs.const_empty_list.copy());

    pallet1[lextype_unknown]        = RGB_COLOR(0x00,0x00,0x00);
    pallet1[lextype_symbol]         = RGB_COLOR(0xff,0xff,0x00);
    pallet1[lextype_symbol_1st]     = RGB_COLOR(0xff,0xff,0x00);
    pallet1[lextype_number]         = RGB_COLOR(0x30,0x18,0x00);
    pallet1[lextype_number_1st]     = RGB_COLOR(0x30,0x18,0x00);
    pallet1[lextype_opinf]          = RGB_COLOR(0x00,0x00,0x00);
    pallet1[lextype_opinf_1st]      = RGB_COLOR(0x00,0x00,0x00);
    pallet1[lextype_oppost]         = RGB_COLOR(0x00,0x00,0x00);
    pallet1[lextype_oppost_1st]     = RGB_COLOR(0x00,0x00,0x00);
    pallet1[lextype_oppre]          = RGB_COLOR(0x00,0x00,0x00);
    pallet1[lextype_oppre_1st]      = RGB_COLOR(0x00,0x00,0x00);
    pallet1[lextype_string]         = RGB_COLOR(0x58,0x60,0x68);
    pallet1[lextype_string_1st]     = RGB_COLOR(0x58,0x60,0x68);
    pallet1[lextype_blank]          = RGB_COLOR(0xd8,0x46,0x00);
    pallet1[lextype_blank_1st]      = RGB_COLOR(0xd8,0x46,0x00);
    pallet1[lextype_pattern_1st]    = RGB_COLOR(0xff,0xff,0x00);
    pallet1[lextype_bracket_open]   = RGB_COLOR(0x00,0x00,0x00);
    pallet1[lextype_bracket_close]  = RGB_COLOR(0x00,0x00,0x00);
    pallet1[lextype_parenth_open]   = RGB_COLOR(0x00,0x00,0x00);
    pallet1[lextype_parenth_close]  = RGB_COLOR(0x00,0x00,0x00);
    pallet1[lextype_curly_open]     = RGB_COLOR(0x00,0x00,0x00);
    pallet1[lextype_curly_close]    = RGB_COLOR(0x00,0x00,0x00);
    pallet1[lextype_expr]           = RGB_COLOR(0x00,0x00,0x00);
    pallet1[lextype_whitespace]     = RGB_COLOR(0x00,0x00,0x00);
    pallet1[lextype_slot]           = RGB_COLOR(0xb0,0x00,0xc0);
    pallet1[lextype_slot_1st]       = RGB_COLOR(0xb0,0x00,0xc0);
    pallet1[lextype_comma]          = RGB_COLOR(0x00,0x00,0x00);
    pallet1[lextype_comment]        = RGB_COLOR(0xa0,0xa8,0xb0);
    pallet1[lextype_comment_1st]    = RGB_COLOR(0xa0,0xa8,0xb0);
    pallet1[lextype_message_name]   = RGB_COLOR(0x58,0x60,0x68);

    pallet1[lextype_MAX + semtype_unknown]  = RGB_COLOR(0x00,0x00,0xc8);
    pallet1[lextype_MAX + semtype_defined]  = RGB_COLOR(0x00,0x00,0x00);
    pallet1[lextype_MAX + semtype_capture]  = RGB_COLOR(0xc0,0x38,0x00);
    pallet1[lextype_MAX + semtype_table]    = RGB_COLOR(0,120,180);
    pallet1[lextype_MAX + semtype_block]    = RGB_COLOR(0x18,0x60,0x40);
    pallet1[lextype_MAX + semtype_module]   = RGB_COLOR(75,0,130);

    cCellForegroundInput        = RGB_COLOR(0x00,0x00,0x00);
    cCellForegroundOutput       = RGB_COLOR(0x00,0x00,0x00);
    cCellForegroundPrint        = RGB_COLOR(0x00,0x00,0x00);
    cCellForegroundMessage      = RGB_COLOR(0xb0,0x50,0x18);
    cCellForegroundText         = RGB_COLOR(0x00,0x00,0x00);
    cCellForegroundSubsubsection= RGB_COLOR(0xB6,0x37,0x08);
    cCellForegroundSubsection   = RGB_COLOR(0xC6,0x6B,0x29);
    cCellForegroundSection      = RGB_COLOR(0xC2,0x4B,0x15);
    cCellForegroundTitle        = RGB_COLOR(0xcc,0x0a,0x02);

    cBackground             = RGB_COLOR(0xff,0xff,0xff);
    cCellBracket            = RGB_COLOR(0x90,0x80,0xc0);
    cCellLabel              = RGB_COLOR(0x30,0x50,0x70);
    cCursor                 = RGB_COLOR(0xfa,0xaa,0x3c);
    cSelectionBackground    = RGB_COLOR(0x33,0x99,0xff);
    cSelectionForeground    = RGB_COLOR(0xff,0xff,0xff);

    cBracketMatch    = RGB_COLOR(0xe0,0xff,0x80);
    cBracketMismatch = RGB_COLOR(0xff,0xd0,0xd0);
};


notebook::~notebook()
{
    delete root;
    if (copy_paste_buffer)
        delete copy_paste_buffer;
};

void notebook::print_cell(cellbox* c)
{
    root->print_cell(c);
    cursor_needs_fitting = true;
}

void notebook::key_insert_char(char16_t c)
{
    cursor_needs_fitting = true;
    root->insert_char(c);
}

void notebook::key_move(moveArg m)
{
    cursor_needs_fitting = true;
    boxbase* b = nullptr;
    root->move(b, m);
    assert(b == nullptr);
}

void notebook::key_insert(insertArg m)
{
    cursor_needs_fitting = true;
    boxbase* b = nullptr;
    root->insert(b, m);
    assert(b == nullptr);
}

void notebook::key_remove(removeArg m)
{
    cursor_needs_fitting = true;
    boxbase* b = nullptr;
    root->remove(b, m);
    assert(b == nullptr);
}

void notebook::key_insert_char16_array(const char16_t * s, size_t len)
{
    cursor_needs_fitting = true;
    boxbase * stuff = new rowbox(s, len);
    root->key_paste(stuff);
    delete stuff;
}

void clipboard_get_data_append(std::string&s);

void clipboard_set_data(const char*s, size_t len);

void notebook::key_paste()
{
    cursor_needs_fitting = true;

    std::string s;
    clipboard_get_data_append(s);
    if (!s.empty() && s != text_buffer)
    {
        text_buffer = s;
        rowbox * stuff = new rowbox(s.c_str(), s.size());
        if (copy_paste_buffer != nullptr)
            delete copy_paste_buffer;
        copy_paste_buffer = stuff;
    }
    root->key_paste(copy_paste_buffer);
}

void notebook::key_copy()
{
    boxbase* b = nullptr;
    root->key_copy(b);
    if (b != nullptr)
    {
        if (copy_paste_buffer != nullptr)
            delete copy_paste_buffer;
        copy_paste_buffer = b;
        b = nullptr;

        wex e(copy_paste_buffer->get_ex());
        text_buffer = ex_tostring_full(e.get());
        clipboard_set_data(text_buffer.data(), text_buffer.size());
    }
}

void notebook::key_toggle_cell_expr()
{
    cursor_needs_fitting = true;
    root->toggle_cell_expr();
}

void notebook::key_makecell(const char * s)
{
    cursor_needs_fitting = true;
    root->makecell(s);
}


void notebook::key_shiftenter()
{
    cursor_needs_fitting = true;
    root->key_shiftenter();
    return;
}

void notebook::mouse_click(int x, int y)
{

}

void notebook::mouse_doubleclick(int x, int y)
{

}

void notebook::zoom_in()
{
    set_zoom128(zoom128 + 1);
}

void notebook::zoom_out()
{
    set_zoom128(zoom128 - 1);
}

void notebook::set_zoom128(int zi)
{
    zi = std::max(zi, 128*MIN_FONT_SIZE/11+1);
    zi = std::min(zi, 128*MAX_FONT_SIZE/11+10);
    zoom128 = zi;
    magnification = 1.0/128*zoom128;
    cursor_needs_fitting = true;
    root->visit(visitarg_InvalidateAll);
}

void notebook::resize(int32_t sizex, int32_t sizey)
{
    if (window_sizex == sizex && window_sizey == sizey)
        return;
    window_sizex = sizex;
    window_sizey = sizey;
    cursor_needs_fitting = true;
    root->visit(visitarg_InvalidateAll);
}

void notebook::get_cursor_window_pos(int32_t & xx, int32_t & yy)
{
    aftransform T;

    root->get_cursor(&T);

    int32_t x = T.orig_x + offx;
    int32_t y = T.orig_y + offy;

    x = std::max(x, (int32_t)0);
    xx = std::min(x, (int32_t)glb_image.sizex - 1);

    y = std::max(y, (int32_t)0);
    yy = std::min(y, (int32_t)glb_image.sizey - 1);
}

void notebook::fitcursor()
{
    if (!cursor_needs_fitting)
        return;

    aftransform T;
    root->get_cursor(&T);

    double x = T.orig_x;
    double y = T.orig_y;

    if (offx < 10 - x)
    {
        offx = 10 - x;
    }
    else if (offx >= glb_image.sizex - x - 10)
    {
        offx = glb_image.sizex - x - 10;
    }
    offx = std::min(offx, 0.0);

    double pady = std::max(20.0, 0.125*glb_image.sizey);

    if (offy < pady - y)
    {
        offy = pady - y;
    }
    else if (offy >= glb_image.sizey - y - pady)
    {
        offy = glb_image.sizey - y - pady;
    }

    cursor_needs_fitting = false;

    _fix_offy();
}

void notebook::measure()
{
    resize(glb_image.sizex, glb_image.sizey);
    default_style.init_const1();
    root->measure(boxmeasurearg(&default_style, glb_image.sizex, 0));
    fitcursor();
}





void notebook::draw_bitmap()
{
//std::cout << "-------------" << std::endl << "measuring notebook:" << std::endl; root->print(0,0,0);

    measure();

//std::cout << "drawing notebook:" << std::endl; root->print(0,offx,offy);

    chardrawcount = 0;
    std::memset(glb_image.pixels, 255, glb_image.pixel_height*glb_image.rowstride);

    aftransform T;

    T.orig_x = 0.0;
    T.orig_y = 0.0;
    T.theta = 0;
    T.cos_theta = 1.0;
    T.sin_theta = 0.0;
//printf("draw_pre\n");
    root->draw_pre(boxdrawarg(this, offx, offy, 0, &T));
//printf("draw_main\n");
    root->draw_main(boxdrawarg(this, offx, offy, 0, &T));
//printf("draw_post\n");
    root->draw_post(boxdrawarg(this, offx, offy, 0, &T));
//printf("draw done\n");

/*
    double q[10][2] = {
        {10.0,10.0},
        {2*300.0,10.0},
        {2*300.0,2*300.0},
        {10.0,2*300.0},
        {10.0,10.0},
    };

    renderPath(&glb_image, q[0], 5, RGB_COLOR(0,162,220));
*/
/*
    for (int i = 0; i<glb_image.pixel_height*glb_image.rowstride; i++)
        glb_image.pixels[i] ^= 255;
*/
#if 0
    svgGlyph _king(64.0, "1.5/2150,0,0,-1.5/2150,1.5/32,36.0/32 "
" M1904 738q-12 -164 -133 -327q-87 -117 -172 -172q-5 -76 -25 -436q-62 -103 -266 -137q-108 -18 -308 -18t-308 18q-204 34 -266 137l-25 436q-122 80 -212 232q-95 158 -95 297q0 136 103 224q97 85 235 85q156 0 307 -112v393h522v-393q151 112 305 112q149 0 247 -94q101 -98 91 -245z");

    svgGlyph wking("1.5/2150,0,0,-1.5/2150,1.5/32,36.0/32 "
"M1190 946v297l-146 -148z"
"M1154 1282h-308l154 -146z"
"M1828 736q7 84 -45 160q-71 104 -219 103q-251 -2 -419 -288q-70 -119 -118 -286q107 -2 205 -15q226 -30 331 -102q247 198 265 428z"
"M1150 911l-150 146l-150 -146h300z"
"M956 1095l-146 148v-297z"
"M1159 840h-318q87 -100 159 -269q74 172 159 269z"
"M973 425q-48 167 -118 286q-168 286 -419 288q-149 1 -220 -103q-52 -77 -44 -160q22 -230 265 -428q157 109 536 117z"
"M1520 243q-222 107 -520 107t-520 -107l9 -129h114l-96 -111v-98q105 108 493 108t493 -108v98l-96 111h114z"
"M1474 -170q0 62 -213 92q-135 19 -261 19t-261 -19q-213 -30 -213 -92t213 -93q138 -20 261 -20t261 20q213 31 213 93z"
"M1117 174q0 -50 -117 -50t-117 50t117 50t117 -50z"
"M1904 738q-12 -164 -133 -327q-87 -117 -172 -172q-5 -76 -25 -436q-62 -103 -266 -137q-108 -18 -308 -18t-308 18q-204 34 -266 137l-25 436q-122 80 -212 232q-95 158 -95 297q0 136 103 224q97 85 235 85q156 0 307 -112v393h522v-393q151 112 305 112q149 0 247 -94q101 -98 91 -245z");

    svgGlyph bking("1.5/2150,0,0,-1.5/2150,1.5/32,36.0/32 "
"M1904 738q-12 -164 -133 -327q-87 -117 -172 -172q-5 -76 -25 -436q-62 -103 -266 -137q-108 -18 -308 -18t-308 18q-204 34 -266 137l-25 436q-122 80 -212 232q-95 158 -95 297q0 136 103 224q97 85 235 85q156 0 307 -112v393h522v-393q151 112 305 112q149 0 247 -94q101 -98 91 -245z"
"M1186 1005q17 20 -2 36l-105 93l102 92q18 16 1 36l-32 37q-18 21 -39 2l-111 -100l-111 100q-21 19 -39 -2l-32 -37q-17 -20 1 -36l102 -92l-105 -93q-19 -16 -2 -36l32 -38q15 -18 37 -4l117 103l117 -103q22 -14 37 4z"
"M1109 763q0 45 -32 74.5t-77 29.5t-77 -29.5t-32 -74.5t32 -74.5t77 -29.5t77 29.5t32 74.5z"
"M1505 55q8 30 -52 46q-56 15 -67 -14q-7 -32 52 -48q58 -16 67 16z"
"M1674 656q3 151 -181 155q-80 -9 -148 -52q-165 -104 -251 -363h-188q-128 385 -399 415q-184 -4 -181 -155q32 -209 207 -325l-51 -42q-5 -12 -5 -26q-1 -34 33 -49q273 79 490 79t490 -79q34 15 33 49q0 14 -5 26l-51 42q174 116 207 325z"
"M1102 145q0 54 -102 54t-102 -54t102 -54t102 54z"
"M1000 -16q-137 0 -266 -14q-220 -24 -220 -76q0 -16 16.5 -31.5t31.5 -11.5q187 53 339 53q50 0 99 -5q49 5 99 5"
"q152 0 339 -53q15 -4 31.5 11.5t16.5 31.5q0 42 -148.5 66t-337.5 24z"
"M614 87q-11 29 -67 14q-60 -16 -52 -46q9 -32 67 -16q59 16 52 48z"
"M1538 623q2 -61 -175 -254q-59 14 -146 27q101 284 274 284q45 0 47 -57z"
"M783 396q-41 -6 -146 -27q-177 193 -175 254q2 63 60 56q53 -7 92 -32q104 -68 169 -251z");

    svgGlyph _queen(64.0, "1.5/2150,0,0,-1.5/2150,1.5/32,36.0/32 "
"M1905 813q0 -45 -31.5 -77.5t-76.5 -32.5l-200 -439l-46 -475q-186 -136 -551 -136q-397 0 -541 136l-56 475l-190 439q-45 0 -76.5 32.5t-31.5 77.5t31.5 77t76.5 32q46 0 77.5 -31.5t31.5 -77.5q0 -30 -15 -54l242 -194l138 524q-88 30 -88 119q0 53 38 89t91 36"
"q55 0 92.5 -35t37.5 -90q0 -47 -34 -84l176 -357l186 357q-34 37 -34 84q0 55 37.5 90t92.5 35q56 0 93.5 -35t37.5 -90q0 -89 -90 -119l128 -524l252 194q-15 24 -15 54q0 46 31.5 77.5t77.5 31.5q45 0 76.5 -32t31.5 -77z");

    svgGlyph wqueen("1.5/2150,0,0,-1.5/2150,1.5/32,36.0/32 "
"M1905 813q0 -45 -31.5 -77.5t-76.5 -32.5l-200 -439l-46 -475q-186 -136 -551 -136q-397 0 -541 136l-56 475l-190 439q-45 0 -76.5 32.5t-31.5 77.5t31.5 77t76.5 32q46 0 77.5 -31.5t31.5 -77.5q0 -30 -15 -54l242 -194l138 524q-88 30 -88 119q0 53 38 89t91 36q55 0 92.5 -35t37.5 -90q0 -47 -34 -84l176 -357l186 357q-34 37 -34 84q0 55 37.5 90t92.5 35q56 0 93.5 -35t37.5 -90q0 -89 -90 -119l128 -524l252 194q-15 24 -15 54q0 46 31.5 77.5t77.5 31.5q45 0 76.5 -32t31.5 -77z"
"M1847 807q0 52 -52 52t-52 -52t52 -52t52 52z"
"M1350 1208q0 65 -68 65q-65 0 -65 -65q0 -67 65 -67q68 0 68 67z"
"M797 1208q0 65 -66 65q-27 0 -47 -19t-20 -46q0 -67 67 -67q66 0 66 67z"
"M1728 699l-308 -242l-157 609l-43 11l-220 -447l-214 445l-37 -9l-168 -609l-314 262l188 -407q268 107 545 107q297 0 545 -107z"
"M1517 244q-218 109 -517 109q-296 0 -516 -109l9 -130h114l-98 -110v-98q163 103 491 98q356 -5 492 -99v99l-98 110h113z"
"M263 813q0 52 -52 52t-52 -52t52 -52t52 52z"
"M1469 -170q0 62 -211 93q-136 20 -258 20q-123 0 -260 -20q-212 -31 -212 -93t212 -92q141 -20 260 -20q118 0 258 20q211 30 211 92z"
"M1113 170q0 -45 -113 -45q-112 0 -112 45t112 45q113 0 113 -45z");

    svgGlyph bqueen("1.5/2150,0,0,-1.5/2150,1.5/32,36.0/32 "
"M1905 813q0 -45 -31.5 -77.5t-76.5 -32.5l-200 -439l-46 -475q-185 -136 -551 -136q-397 0 -541 136l-56 475l-190 439q-45 0 -76.5 32.5t-31.5 77.5t31.5 77t76.5 32q46 0 77.5 -31.5t31.5 -77.5q0 -30 -15 -54l242 -194l138 524q-88 29 -88 119q0 53 38 89t91 36q55 0 92.5 -35t37.5 -90q0 -47 -34 -84l176 -357l186 357q-34 37 -34 84q0 55 37.5 90t92.5 35q56 0 93.5 -35t37.5 -90q0 -89 -90 -119l128 -524l252 194q-15 24 -15 54q0 46 31.5 77.5t77.5 31.5q45 0 76.5 -32t31.5 -77z"
"M1506 45q0 26 -53 40q-68 17 -68 -23q0 -25 53 -40q68 -20 68 23z"
"M1513 275q-14 29 -131 64q-151 45 -382 45t-382 -45q-117 -35 -131 -64q-8 -16 -8 -32q0 -41 33 -47q-9 2 133 39q154 40 355 40t355 -40q142 -37 133 -39q33 6 33 47q0 16 -8 32z"
"M1104 133q0 54 -104 54t-104 -54q0 -56 104 -56t104 56z"
"M1000 -9q-140 0 -249 -15q-203 -28 -203 -97q0 -42 30 -47q-9 2 100 29q122 30 322 30t322 -30q109 -27 100 -29q30 5 30 47q0 69 -203 97q-109 15 -249 15zM553 85q-53 -14 -53 -40q0 -43 68 -23q53 15 53 40q0 40 -68 23z");

    svgGlyph _rook(56.0, "1.5/2150,0,0,-1.5/2150,1.5/32,36.0/32 "
"M399 992v333h234v-154h210v155h314v-155h210v154h234v-333l-156 -155v-659l129 -132v-52l106 -100v-226h-1360v226l106 100v52l129 132v659z");

    svgGlyph wrook("1.5/2150,0,0,-1.5/2150,1.5/32,36.0/32 "
"M399 992v333h234v-154h210v155h314v-155h210v154h234v-333l-156 -155v-659l129 -132v-52l106 -100v-226h-1360v226l106 100v52l129 132v659z"
"M635 211h730v603h-730v-603z"
"M476 1013l130 -124h784l134 124v239h-82v-150h-363v150h-158v-150h-363v150h-82v-239z"
"M510 46v-52l-107 -100v-148h1196v148l-109 100v52l-96 94h-788z");

    svgGlyph brook("1.5/2150,0,0,-1.5/2150,1.5/32,36.0/32 "
"M1680 -332h-1360v226l106 100v52l129 132v659l-156 155v333h234v-154h210v155h314v-155h210v154h234v-333l-156 -155v-659l129 -132v-52l106 -100v-226z"
"M1363 780v133h-730v-133h730z"
"M1363 106v133h-730v-133h730z");

    svgGlyph _wknight(64.0, "1.5/2150,0,0,-1.5/2150,1.5/32,36.0/32 "
"M1845 66q-6 -333 -56 -396l-1161 -2q-64 73 -64 152q0 107 90 178q313 248 306 325q-87 -16 -160 -16q-53 0 -98 8q-102 -205 -236 -207q-12 -3 -23 -3q-36 0 -55 30q-23 -31 -56 -31q-7 0 -17 2q-85 17 -121 54q-40 42 -40 125q0 16 1 30q7 101 124 222q123 127 135 195q11 60 19 76q24 46 112 127q-5 8 -18 113q-15 121 -15 196q0 23 2 41q87 -66 125 -100q7 141 73 221q39 -78 150 -183q372 8 623 -250q192 -197 291 -518q72 -232 69 -389z");

    svgGlyph wknight("1.5/2150,0,0,-1.5/2150,1.5/32,36.0/32 "
"M1845 66q-6 -333 -56 -396l-1161 -2q-64 73 -64 152q0 107 90 178q313 248 306 325q-87 -16 -160 -16q-53 0 -98 8q-102 -205 -236 -207q-12 -3 -23 -3q-36 0 -55 30q-23 -31 -56 -31q-7 0 -17 2q-85 17 -121 54q-40 42 -40 125q0 16 1 30q7 101 124 222q123 127 135 195q11 60 19 76q24 46 112 127q-5 8 -18 113q-15 121 -15 196q0 23 2 41q87 -66 125 -100q7 141 73 221q39 -78 150 -183q372 8 623 -250q192 -197 291 -518q72 -232 69 -389z"
"M652 1101q-14 50 -51 40q-41 -11 -27 -62q13 -49 54 -38q38 10 24 60z"
"M1766 154l-94 98l73 125l-140 94l52 130l-133 85l29 117l-121 19l-21 113l-117 16l-27 88l-75 -15l-19 63l-110 -11l-36 71l-138 -29q-68 38 -91 57q-10 8 -73 71q-22 -46 -22 -84q0 -42 46 -95q32 -37 32 -74q0 -16 -6 -31l-107 73q-21 -21 -79 -132q-77 -72 -108 -135q-22 -45 -44 -138q-24 -58 -184 -213q-42 -41 -42 -127q0 -103 113 -107q84 101 130 101q5 0 8 -1q-18 -32 -18 -62q0 -18 7 -32q11 -5 26 -5q83 0 168 139q-43 41 -43 86q0 7 1 14q70 -36 218 -36q207 0 242 82q-14 -158 -88 -269q-82 -123 -307 -296q-40 -31 -40 -90q0 -39 26 -69q262 30 544 30q284 0 510 -34l52 188l-80 85z"
"M848 732q-71 -36 -146 -73q-82 -8 -84 27q-3 46 19.5 71.5t68.5 22.5q56 -3 142 -48z"
"M399 354q9 -23 -40 -41q-45 -17 -56 6q-7 23 40 40q47 18 56 -5z");

    svgGlyph _bknight(64.0, "1.5/2150,0,0,-1.5/2150,1.5/32,36.0/32 "
"M1767 292q13 -81 13 -175q0 -212 -57 -437h-1062q-75 0 -75 58q0 50 52 131q43 67 98 122q277 118 277 451q0 19 -1 38q-154 -169 -242 -169q-19 0 -34 10q-104 -191 -192 -228q-60 -25 -89 -25q-25 0 -25 18q0 21 20 37q6 5 88 59q38 25 46 62q4 19 4 34q0 46 -38 46q-54 0 -92 -47q-34 -41 -41 -98q-45 -68 -85 -68q-30 0 -89 73q-42 52 -42 115q0 66 32 97q75 73 129 155q2 3 76 123q5 173 146 338q-58 104 -58 229q0 39 6 69l106 -106h40l10 180q79 -119 122 -155q43 -26 85 -52q193 -39 229 -51q538 -182 643 -834z");

    svgGlyph bknight("1.5/2150,0,0,-1.5/2150,1.5/32,36.0/32 "
"M1767 292q13 -81 13 -175q0 -212 -57 -437h-1062q-75 0 -75 58q0 50 52 131q43 67 98 122q277 118 277 451q0 19 -1 38q-154 -169 -242 -169q-19 0 -34 10q-104 -191 -192 -228q-60 -25 -89 -25q-25 0 -25 18q0 21 20 37q6 5 88 59q38 25 46 62q4 19 4 34q0 46 -38 46q-54 0 -92 -47q-34 -41 -41 -98q-45 -68 -85 -68q-30 0 -89 73q-42 52 -42 115q0 66 32 97q75 73 129 155q2 3 76 123q5 173 146 338q-58 104 -58 229q0 39 6 69l106 -106h40l10 180q79 -119 122 -155q43 -26 85 -52q193 -39 229 -51q538 -182 643 -834z"
"M830 995q-46 10 -57 29q-14 39 -27 78q-47 -63 -47 -124q0 -29 8 -64q11 -48 70 -48q13 0 30 4q34 8 34 50q0 28 -11 75z"
"M1683 31q0 278 -48 445q-55 191 -245 389q-201 209 -363 212q-41 1 -41 -43q0 -27 32 -39q124 -47 200 -109q38 -31 159 -160q108 -115 162 -320q45 -170 45 -361q0 6 -8 -152q-6 -118 5 -135h69q12 16 23 127q10 100 10 146z"
"M746 630q1 5 1 12q0 56 -40.5 90t-97.5 34q-30 0 -61 -7q-62 -14 -62 -51q0 -11 14 -17l54 -22q19 -13 23 -62q68 31 109 31q22 0 60 -8z"
"M444 340q4 16 -30 28q-28 10 -47 10q-93 0 -93 -101q0 -50 40 -50q20 0 32 17q27 38 63 52q26 10 35 44z");

    svgGlyph _pawn(54.0, "1.5/2150,0,0,-1.5/2150,1.5/32,36.0/32 M1600 -295h-1200q-5 157 122 228q99 55 268 54l6 85q-62 34 -62 135q0 127 108 326v131h-361q-8 59 31 134q62 120 220 163q69 42 98 52q12 4 25 19q-21 32 -26 43q-13 27 -14 47q-2 74 54 125.5t131 51.5q74 0 130.5 -51.5t54.5 -125.5q-2 -67 -40 -90q10 -41 123 -71"
"q254 -68 254 -266q0 -20 -3 -31h-361v-131q109 -199 109 -326q0 -101 -63 -135l6 -85q230 2 329 -101q63 -66 63 -157q0 -15 -2 -24z");

    svgGlyph wpawn("1.5/2150,0,0,-1.5/2150,1.5/32,36.0/32 "
"M1600 -295h-1200q-5 157 122 228q99 55 268 54l6 85q-62 34 -62 135q0 127 108 326v131h-361q-8 59 31 134q62 120 220 163q69 42 98 52q12 4 25 19q-21 32 -26 43q-13 27 -14 47q-2 74 54 125.5t131 51.5q74 0 130.5 -51.5t54.5 -125.5q-2 -67 -40 -90q10 -41 123 -71q254 -68 254 -266q0 -20 -3 -31h-361v-131q109 -199 109 -326q0 -101 -63 -135l6 -85q230 2 329 -101q63 -66 63 -157q0 -15 -2 -24z"
"M1423 732q12 108 -189 167q-124 37 -185 129q52 29 52 94q0 43 -29 71.5t-72 28.5t-72.5 -28.5t-28.5 -71.5q1 -75 52 -94q-57 -95 -185 -129q-70 -19 -122 -57q-68 -50 -67 -110h846z"
"M1494 -218q-14 50 -47 75q-83 63 -294 65l-26 163q26 7 36 21q16 24 26 95q-2 48 -13 88q-27 99 -100 232v143h-152v-143q-124 -219 -113 -320q7 -63 30 -94q5 -7 32 -22l-26 -163q-213 -2 -299 -65q-50 -37 -46 -75h992z");

    svgGlyph bpawn("1.5/2150,0,0,-1.5/2150,1.5/32,36.0/32 "
"M1600 -295h-1200q-5 157 122 228q99 55 268 54l6 85q-62 34 -62 135q0 127 108 326v131h-361q-8 59 31 134q62 120 220 163q69 42 98 52q12 4 25 19q-21 32 -26 43q-13 27 -14 47q-2 74 54 125.5t131 51.5q74 0 130.5 -51.5t54.5 -125.5q-2 -67 -40 -90q10 -41 123 -71q254 -68 254 -266q0 -20 -3 -31h-361v-131q109 -199 109 -326q0 -101 -63 -135l6 -85q230 2 329 -101q63 -66 63 -157q0 -15 -2 -24z");


    svgGlyph _wbishop(64.0, "1.5/2150,0,0,-1.5/2150,1.5/32,36.0/32 "
"M1902 -401q-126 60 -393 60q-121 0 -294 -8q-49 107 -213 107q-168 0 -217 -107q-173 8 -294 8q-267 0 -393 -60q-17 36 -17 90q0 217 315 217q103 0 237 -25q-23 96 -23 183q0 43 6 84q-41 59 -78 149q-57 138 -57 261q0 241 312 608q-50 24 -50 96q0 46 32.5 78.5t78.5 32.5q47 0 79.5 -32t33.5 -79q1 -64 -57 -96q47 -59 94 -94q78 65 96 107q-40 40 -40 83q0 47 33.5 79t80.5 32q49 0 79.5 -31.5t32.5 -80.5q3 -61 -48 -82q253 -307 250 -597q-1 -128 -38 -264q-33 -122 -66 -170q6 -41 6 -84q0 -87 -23 -183q134 25 237 25q315 0 315 -217q0 -54 -17 -90z");

    svgGlyph wbishop("1.5/2150,0,0,-1.5/2150,1.5/32,36.0/32 "
"M1902 -401q-125 59 -393 59q-121 0 -294 -7q-48 107 -213 107q-169 0 -217 -107q-173 7 -294 7q-268 0 -393 -59q-18 36 -18 90q0 218 317 218q102 0 236 -26q-23 97 -23 183q0 43 6 84q-42 59 -79 148q-57 138 -57 261q0 8 2 35q17 227 311 574q-50 23 -50 96q0 46 32.5 78.5t78.5 32.5q47 0 79.5 -32t33.5 -79q1 -64 -57 -96q46 -60 94 -94q78 65 96 107q-40 40 -40 83q0 47 33.5 79t80.5 32q48 0 79 -32t33 -80q3 -60 -48 -82q254 -307 250 -597q-2 -122 -38 -263q-32 -124 -66 -171q6 -41 6 -84q0 -86 -23 -183q134 26 237 26q316 0 316 -218q0 -54 -18 -90z"
"M1216 1267q0 45 -42 45q-44 0 -44 -45q0 -44 44 -44q42 0 42 44z"
"M896 1267q0 45 -42 45q-44 0 -44 -45q0 -44 44 -44q42 0 42 44z"
"M1413 569q2 148 -70 303q-51 110 -169 274q-54 -61 -132 -140q37 -40 84 -94q226 -259 245 -564q40 72 42 221z"
"M1309 394q-8 112 -101 279q-124 223 -354 447q-102 -71 -203 -269t-101 -324q0 -150 118 -279q96 37 182 47q43 5 150 5q177 0 280 -56q36 54 29 150z"
"M1856 -297q1 8 1 16q0 119 -250 119q-180 0 -453 -60q34 -52 354 -63z"
"M1309 -42l-4 79l-53 46l57 34v50q-73 58 -309 58t-309 -58v-50l57 -34l-53 -46l-4 -79q106 31 309 31t309 -31zM1265 -123q0 52 -265 52t-265 -52q0 -51 265 -51t265 51zM846 -222q-273 60 -453 60q-250 0 -250 -119q0 -8 1 -16l348 12q320 11 354 63zM1069 102q0 -36 -69 -36t-69 36q0 35 69 35t69 -35z");

    svgGlyph _bbishop(64.0, "1.5/2150,0,0,-1.5/2150,1.5/32,36.0/32 "
"M1902 -401q-126 60 -393 60q-121 0 -294 -8q-49 107 -213 107q-168 0 -217 -107q-173 8 -294 8q-267 0 -393 -60q-17 36 -17 90q0 217 315 217q103 0 237 -25q-23 96 -23 183q0 43 6 84q-41 59 -78 149q-57 138 -57 261q0 241 312 608q-50 24 -50 96q0 46 32.5 78.5t78.5 32.5q47 0 79.5 -32t33.5 -79q1 -64 -57 -96q47 -59 94 -94q78 65 96 107q-40 40 -40 83q0 47 33.5 79t80.5 32q49 0 79.5 -31.5t32.5 -80.5q3 -61 -48 -82q253 -307 250 -597q-1 -128 -38 -264q-33 -122 -66 -170q6 -41 6 -84q0 -87 -23 -183q134 25 237 25q315 0 315 -217q0 -54 -17 -90z");

    svgGlyph bbishop("1.5/2150,0,0,-1.5/2150,1.5/32,36.0/32 "
"M1902 -401q-126 60 -393 60q-121 0 -294 -8q-49 107 -213 107q-168 0 -217 -107q-173 8 -294 8q-267 0 -393 -60q-17 36 -17 90q0 217 315 217q103 0 237 -25q-23 96 -23 183q0 43 6 84q-41 59 -78 149q-57 138 -57 261q0 241 312 608q-50 24 -50 96q0 46 32.5 78.5t78.5 32.5q47 0 79.5 -32t33.5 -79q1 -64 -57 -96q47 -59 94 -94q78 65 96 107q-40 40 -40 83q0 47 33.5 79t80.5 32q49 0 79.5 -31.5t32.5 -80.5q3 -61 -48 -82q253 -307 250 -597q-1 -128 -38 -264q-33 -122 -66 -170q6 -41 6 -84q0 -87 -23 -183q134 25 237 25q315 0 315 -217q0 -54 -17 -90z"
"M1359 470q0 213 -133 421q-62 97 -110 119q-30 14 -56 14q-46 0 -48 -38q-1 -16 69 -113q88 -122 110 -189q17 -61 51 -184q41 -143 75 -143q42 0 42 113z"
"M1309 101q7 18 -27 27q-31 8 -38 -10q-3 -20 29 -29q33 -9 36 12z"
"M1311 266q-46 36 -86 50q-85 30 -225 30q-209 0 -311 -80q-19 -15 -19 -41q0 -50 34 -50l8 2q73 30 97 36q75 20 191 26q155 -8 288 -62q5 -2 8 -2q34 0 34 50q0 26 -19 41z"
"M1082 147q0 30 -82 30t-82 -30q0 -31 82 -31t82 31z"
"M1281 -24q-63 34 -89 43q-79 27 -192 27q-155 0 -281 -70q-22 -12 -22 -39q0 -52 36 -52q5 0 9 2q73 33 94 40q64 22 164 22t164 -22q21 -7 94 -40q4 -2 9 -2q36 0 36 52q0 27 -22 39z"
"M756 118q-7 18 -38 10q-34 -9 -27 -27q3 -21 36 -12q32 9 29 29z");

    double scale = 46;
    double off = 20;
    double epsilon = 0.5/scale;

    char board[65] =
/*
                     "....r..k"
                     ".NR..p.p"
                     ".....q.."
                     "...Q...."
                     "...p...."
                     "PP.n.pPP"
                     ".....P.."
                     "......K.";
*/
                     "rnbqkbnr"
                     "pppppppp"
                     "........"
                     "........"
                     "........"
                     "........"
                     "PPPPPPPP"
                     "RNBQKBNR";


    for (int i = 0; i <= 8; i++)
    {
        renderLine(&glb_image, 0*scale + off, i*scale + off, 8*scale + off, i*scale + off,  2.0, RGB_COLOR(175,175,124));
        renderLine(&glb_image, i*scale + off, 0*scale + off, i*scale + off, 8*scale + off,  2.0, RGB_COLOR(175,175,124));
    }

    for (int i = 0; i < 8; i++)
    for (int j = 0; j < 8; j++)
    {
        affTran W(scale,0,0,scale,scale*j + off,scale*i + off);
        uint32_t color = ((i+j)&1) ? RGB_COLOR(156,156,100) : RGB_COLOR(206,206,156);
//        renderRect(&glb_image, color, W, epsilon, 1-epsilon, epsilon, 1-epsilon);

        switch (board[8*i+j])
        {
        case 'k':
            renderGlyph(&glb_image, 0xFFFFFF, _king, W);
            renderGlyph(&glb_image, 0x000000, bking, W);
            break;
        case 'K':
            renderGlyph(&glb_image, 0xFFFFFF, _king, W);
            renderGlyph(&glb_image, 0x000000, wking, W);
            break;
        case 'q':
            renderGlyph(&glb_image, 0xFFFFFF, _queen, W);
            renderGlyph(&glb_image, 0x000000, bqueen, W);
            break;
        case 'Q':
            renderGlyph(&glb_image, 0xFFFFFF, _queen, W);
            renderGlyph(&glb_image, 0x000000, wqueen, W);
            break;
        case 'r':
            renderGlyph(&glb_image, 0xFFFFFF, _rook, W);
            renderGlyph(&glb_image, 0x000000, brook, W);            
            break;
        case 'R':
            renderGlyph(&glb_image, 0xFFFFFF, _rook, W);
            renderGlyph(&glb_image, 0x000000, wrook, W);            
            break;
        case 'b':
            renderGlyph(&glb_image, 0xFFFFFF, _bbishop, W);
            renderGlyph(&glb_image, 0x000000, bbishop, W);            
            break;
        case 'B':
            renderGlyph(&glb_image, 0xFFFFFF, _wbishop, W);
            renderGlyph(&glb_image, 0x000000, wbishop, W);            
            break;
        case 'n':
            renderGlyph(&glb_image, 0xFFFFFF, _bknight, W);
            renderGlyph(&glb_image, 0x000000, bknight, W);            
            break;
        case 'N':
            renderGlyph(&glb_image, 0xFFFFFF, _wknight, W);
            renderGlyph(&glb_image, 0x000000, wknight, W);            
            break;
        case 'p':
            renderGlyph(&glb_image, 0xFFFFFF, _pawn, W);
            renderGlyph(&glb_image, 0x000000, bpawn, W);            
            break;
        case 'P':
            renderGlyph(&glb_image, 0xFFFFFF, _pawn, W);
            renderGlyph(&glb_image, 0x000000, wpawn, W);            
            break;
        }
    }
#endif
}


void notebook::_fix_offy()
{
    offy = std::max(offy, -root->sizey + (double)glb_image.sizey/2);
    offy = std::min(offy, 0.0);
}

void notebook::scroll_up()
{
    offy += glb_image.sizey/10;
    _fix_offy();
}

void notebook::scroll_down()
{
    offy -= glb_image.sizey/10;
    _fix_offy();
}


void notebook::key_pageup()
{
    offy += glb_image.sizey;
    _fix_offy();
}

void notebook::key_pagedown()
{
    offy -= glb_image.sizey;
    _fix_offy();
}

void notebook::save(const char * fstr)
{
    FILE *fp = fopen(fstr,"wb");
    if (!fp)
    {
        std::cerr << "save failed" << std::endl;
        return;
    }
    wex e(root->get_ex());
    std::string s = ex_tostring_full(e.get());
    fwrite(s.c_str(), 1, s.size(), fp);
    fclose(fp);
    if (filestring.compare(fstr) != 0)
    {
        filestring.assign(fstr);
    }
}

bool notebook::open(const char * fstr)
{
    std::ifstream is;
    is.open(fstr, std::ios::in);
    if ((is.rdstate() & std::ifstream::failbit) != 0)
    {
        std::cerr << "failed to open notebook" << fstr << std::endl;
        return false;
    }

    std::vector<uex> v;
    syntax_report sr;
    ex_parse_file(v, is, false, sr);
    if (sr.have_error() || v.size() != 1)
    {
        std::cerr << "syntax error in notebook " << fstr << std::endl;
        return false;
    }

    boxbase* r = boxbase_from_ex(v[0].get());
    if (r->get_type() != BNTYPE_ROOT)
    {
        delete r;
        std::cerr << "improper format of notebook " << fstr << std::endl;
        return false;    
    }

    if (root != nullptr)
        delete root;
    root = dynamic_cast<rootbox*>(r);
    root->nb = this;

    filestring.assign(fstr);
    return true;
}
