#include "globalstate.h"
#include "code.h"
#include "dirichlet.h"
#include "acb_dirichlet.h"
#include "xacb_t.h"
#include "flintarb_wrappers.h"

ex scode_sDirichletCharacter(er e)
{
//std::cout << "scode_sDirichletCharacter: " << e << std::endl;
    assert(eis_node(e));

    if (elength(e) != 1 || !eis_int(echild(e,1)))
        return ecopy(e);

    er e1 = echild(e,1);
    er h = echild(e,0);
    if (!ehas_head_sym_length(h, gs.sym_sDirichletCharacter.get(), 2))
        return ecopy(e);

    er h1 = echild(h,1);
    er h2 = echild(h,2);
    if (!eis_intsm(h1) || !eis_int(h2))
        return ecopy(e);

    ulong q = eintsm_get(h1);
    if (eintsm_get(h1) < 1)
        return ecopy(e);

    ulong m = fmpz_fdiv_ui(eint_data(h2), q);
    if (n_gcd(q, m) != 1)
        return ecopy(e);

    ulong n = fmpz_fdiv_ui(eint_data(e1), q);
    if (n_gcd(q, n) != 1)
        return emake_cint(0);

    dirichlet_group_t G;
    dirichlet_group_init(G, q);
    ulong d = G->expo;

    dirichlet_char_t chi, nchi;
    dirichlet_char_init(chi, G);
    dirichlet_char_log(chi, G, m);
    dirichlet_char_init(nchi, G);
    dirichlet_char_log(nchi, G, n);

    ulong a = 2*dirichlet_pairing_char(G, chi, nchi);
    ulong g = n_gcd(a, d);

    dirichlet_char_clear(chi);
    dirichlet_char_clear(nchi);
    dirichlet_group_clear(G);

    return ex_powx(emake_cint(-1), emake_rat(slong(a/g), slong(d/g)));
}

ex dcode_sDirichletCharacter(er e)
{
//std::cout << "dcode_sDirichletCharacter: " << e << std::endl;
    assert(ehas_head_sym(e, gs.sym_sDirichletCharacter.get()));

    if (enode_statusflags(e) & ESTATUSFLAG_VALIDQ)
        return ecopy(e);

    if (elength(e) != 2)
        return _handle_message_argx(e, 2);

    er e1 = echild(e,1);
    er e2 = echild(e,2);
    if (eis_int(e1))
    {
        if (fmpz_sgn(eint_data(e1)) <= 0)
        {
            _gen_message(echild(e,0), "badq", "Modulus `1` is not positive.", ecopy(e1));
        }
        else if (eis_int(e2))
        {
            if (!fmpz_coprime(eint_data(e1), eint_data(e2)))
                _gen_message(echild(e,0), "bade", "Character exponent `1` is not coprime to modulus `2`.", ecopy(e2), ecopy(e1));

            if (fmpz_sgn(eint_data(e2)) < 0 ||
                fmpz_cmp(eint_data(e2), eint_data(e1)) >= 0)
            {
                fmpz_t ne2;
                fmpz_init(ne2);
                fmpz_fdiv_r(ne2, eint_data(e2), eint_data(e1));
                ex t = emake_int_clear(ne2);
                t = emake_node(ecopychild(e,0), ecopychild(e,1), t);
                enode_statusflags(t) |= ESTATUSFLAG_VALIDQ;
                return t;
            }

            enode_statusflags(e) |= ESTATUSFLAG_VALIDQ;
        }
    }

    return ecopy(e);
}


ex dcode_sDirichletL(er e)
{
//std::cout << "dcode_sDirichletL: " << e << std::endl;
    assert(ehas_head_sym(e, gs.sym_sDirichletL.get()));

    if (elength(e) != 2)
        return _handle_message_argx2(e);

    er c = echild(e,1);
    if (!ehas_head_sym_length(c, gs.sym_sDirichletCharacter.get(), 2))
        return ecopy(e);

    er c1 = echild(c,1);
    er c2 = echild(c,2);
    if (!eis_intsm(c1) || !eis_int(c2))
        return ecopy(e);

    ulong q = eintsm_get(c1);
    if (eintsm_get(c1) < 1)
        return ecopy(e);

    ulong m = fmpz_fdiv_ui(eint_data(c2), q);
    if (n_gcd(q, m) != 1)
        return ecopy(e);

    er S = echild(e,2);
    if (!eis_real(S))
        return ecopy(e);


    dirichlet_group_t G;
    dirichlet_group_init(G, q);
    dirichlet_char_t chi;
    dirichlet_char_init(chi, G);
    dirichlet_char_log(chi, G, m);

    xacb_t res, s;
    arb_set(acb_realref(s.data), ereal_data(S));
    acb_dirichlet_l(res.data, s.data, G, chi, s.wprec() + EXTRA_PRECISION_BASIC);

    dirichlet_char_clear(chi);
    dirichlet_group_clear(G);


    return emake_cmplx_move(res.data);
}
