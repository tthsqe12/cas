

rfq_nmod_t::mul(const ulong* b, const ulong* c)
{
    compute 


}

rfq_nmod_t::reduce(ulong* a)
{
    for (slong j = 0; j < d - 1; j++)
        reduce3(r[j], p[3*(d+j)+2], p[3*(d+j)+1], p[3*(d+j)+0]);

    for (slong i = 0; i < d; i++)
    for (slong j = 0; j < d - 1; j++)
        addmul(p[3*i+2], p[3*i+1], p[3*i+0], r[j], m[d*i+j]);

    for (slong i = 0; i < d; i++)
        reduce3(a[i], p[3*i+2], p[3*i+1], p[3*i+0]);
}