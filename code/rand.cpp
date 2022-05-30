#include "rand.h"
#include "flintarb_wrappers.h"

/*
SeedRandom[
  12360270945421506099765375950419569273198467234936204812865465999525\
7623109708538877860925408394368448624071031081828236607796317969385775\
4897912671510956691707029835002042778409991739027488028285620627477847\
8632228462952969559809870833719448767087136300339834806286644888228761\
0571899085123016173685313316447200603056620208228194860100439722138310\
4293154655188863631049365314881119974826918372516362326115079224728641\
2602221060167318282129317754113704318615849220954166200159997377708233\
7866445015489302515067652349904823849403911289398949202423721198859910\
2767364323934868959398224702943405965315096871694077537921514087681610\
5253976761559026369032969869195427625283948490477079451839108429318466\
9208110362484699255486531915469714968733103778022414453492790273898850\
1926247307827789407032810194538591296856729000279231659368769933023758\
0669013789084144764119524361168735148789766477129650239166728167336673\
2307440421198524141752079224272952354893589517553639391855579155660612\
3065460314530115846021541651330724245166597470485287234352938165171357\
6566170006213986654194869290945648251352923806821727776372235474194867\
2675702580951021977344588193834407079948264958695169763042796760530952\
3591875908898183932777973833537475576963589446045819076661515220069820\
6554058976145500119655852914599606023625634708716713794291555774992670\
5530780983485860086289515624817795242127161164406899938333241509741146\
3687346922684677272848247956182887052568278883670898112653450791441501\
1725312735087031410910243755821253505009053550174393404549699019823569\
3050];
Partition[Table[RandomInteger[], {10000}], 80] // Grid
*/

extca_rand::extca_rand()
{
    startover();
    shuffle();
    off = 0;
    bitsleft = 64;
}

void extca_rand::startover()
{
    D[ 0]=0xa888f71022e181bdULL;D[ 1]=0x2385a85e36087989ULL;D[2]=0x20ff49270ec954cfULL;D[3]=0xa68f432825ff749bULL;
    D[ 4]=0x87bfabfaee5c4d9eULL;D[ 5]=0xfc0d678d213f99cULL;D[6]=0x439271a5c7c22777ULL;D[7]=0xabe5acb9642a8103ULL;
    D[ 8]=0x5c17a5228613598dULL;D[ 9]=0x89de293826b2a148ULL;D[10]=0x4454535aec23e7c6ULL;D[11]=0x7226cb254a47cbc6ULL;
    D[12]=0xba8b66c3835aad7eULL;D[13]=0xfb59214b81d4a98cULL;D[14]=0xe62dc0746ad746beULL;D[15]=0xd1f4194ecb6cd9d5ULL;
    D[16]=0x78701f81f42c1948ULL;D[17]=0x8fccc9aaea1b68e9ULL;D[18]=0x464adfa0a258c07dULL;D[19]=0xf1f3ed907b3c47dULL;
    D[20]=0x8dec694af72add42ULL;D[21]=0xbfee489c7e8ef175ULL;D[22]=0xbfa60a9c4435876eULL;D[23]=0x452d611c7f055daULL;
    D[24]=0x8c8f3854e3fa1448ULL;D[25]=0x9e83bfb8c0b0e837ULL;D[26]=0xa942abdb2d91d070ULL;D[27]=0x6f52de0bf1e47693ULL;
    D[28]=0x34c4d38ad1212667ULL;D[29]=0x4b06768ef21391bdULL;D[30]=0x8e1412e997c9a8fdULL;D[31]=0x50bfda09377ea28aULL;
    D[32]=0x078800b5470871eaULL;D[33]=0xa9aaa916da96c5fbULL;D[34]=0x22152498eca9c50aULL;D[35]=0x6c51da24b87250b0ULL;
    D[36]=0x1d46e78f65a2ef55ULL;D[37]=0x2da0d6c791de8be7ULL;D[38]=0x90295517bc846686ULL;D[39]=0x998a54e3218ac814ULL;
    D[40]=0xd3ec5836257b96f8ULL;D[41]=0x19908f7cac6a00f5ULL;D[42]=0xfaa481ec4c16019cULL;D[43]=0x97e6398dd5dee075ULL;
    D[44]=0x48475ccfba79d982ULL;D[45]=0xd6a1062d4f3bcd4cULL;D[46]=0x87a0959ece65f509ULL;D[47]=0x6719e57ed6ede60dULL;
    D[48]=0xf648daa575d9c6f7ULL;D[49]=0xe5bb2da9c118e114ULL;D[50]=0xaaebc6474fc5a834ULL;D[51]=0x874dcaf472a76dfcULL;
    D[52]=0x56309d44ddfa7109ULL;D[53]=0x9001cac5d3254553ULL;D[54]=0xa46727129321f695ULL;D[55]=0x357e1feae487b4b8ULL;
    D[56]=0xdf6d1e1461a5876fULL;D[57]=0x36fe3f03e1620147ULL;D[58]=0x3bb292e092dbec1eULL;D[59]=0x4b054997f76a743ULL;
    D[60]=0xd4629048376cb213ULL;D[61]=0x06f9bcb003bf6f82ULL;D[62]=0x5d36b34645d974d9ULL;D[63]=0x91a12bd764896c61ULL;
    D[64]=0x6646abd7f82ff67fULL;D[65]=0x90447798b42b05c7ULL;D[66]=0xbc10e5e87050f778ULL;D[67]=0x9ab994748afd5155ULL;
    D[68]=0x05fe8ae3c32bb299ULL;D[69]=0xdfd34d860ea09ccaULL;D[70]=0xb99a9896a7dc966dULL;D[71]=0x9f12f07a7dbb9ec7ULL;
    D[72]=0xe1a3bb6258d768e0ULL;D[73]=0x8a580cc45515cf7dULL;D[74]=0x160c96cc71849d40ULL;D[75]=0x66c6624197f3d708ULL;
    D[76]=0xb1acb50934546050ULL;D[77]=0xae0236b997921d2aULL;D[78]=0xf632a828b5ef628eULL;D[79]=0x14e8210276047ffaULL;
}

void extca_rand::shuffle()
{
    uint64_t E[132];

    for (int i = 0; i < 80; i++)
    {
        E[i] = D[i];
    }
    for (int i = 0; i < 52; i++)
    {
        E[80+i] = 2*D[i] + (D[i]>>63);
    }
    for (int i = 0; i < 80; i++)
    {
        D[i] = ~((E[i+52]|E[i+12]|E[i+29])^E[i]^E[i+5]);
    }
}

void extca_rand::print()
{
    for (int i = 0; i < 80; i += 4)
    {
        printf("0x%016llx 0x%016llx 0x%016llx 0x%016llx\n", D[i+0], D[i+1], D[i+2], D[i+3]);
    }
}

void extca_rand::seed(const unsigned char * data, size_t len)
{
    startover();
    off = 0;
    bitsleft = 64;

    for (size_t i = 0; i < len; i++)
    {
        Db[i%640] ^= data[i];
    }

    for (int i = (len == 0 ? 1 : 197); i > 0; i--)
    {
        shuffle();
    }
}

void extca_rand::getbits(xfmpz_t &x, ulong count)
{
    if (count <= FLINT_BITS)
    {
        fmpz_set_ui(x.data, getbits(count));
        return;
    }

    fmpz_zero(x.data);

    while (count != 0)
    {
        count--;

        if (bitsleft == 0)
        {
            off += 4;
            if (off >= 80)
            {
                shuffle();
                off = 0;
            }
            bitsleft = 64;
        }
        bitsleft -= 1;

        uint64_t m = D[off];
        uint64_t n = m << 1;
        D[off] = n;
        if (n < m)
        {
            fmpz_setbit(x.data, count);
        }
    }
}

ulong extca_rand::getbits(unsigned int count)
{
    assert(0 < count && count <= FLINT_BITS && FLINT_BITS <= 64);

    if (bitsleft == 0)
    {
        off += 4;
        if (off >= 80)
        {
            shuffle();
            off = 0;
        }
        bitsleft = 64;
    }

    if (likely(bitsleft >= count))
    {
        bitsleft -= count;
        ulong r = SHLD(0, D[off], count);
        D[off] = D[off] << count;
        return r;
    }
    else
    {
        count -= bitsleft;
        ulong r = SHLD(0, D[off], bitsleft);

        off += 4;
        if (unlikely(off >= 80))
        {
            shuffle();
            off = 0;
        }

        bitsleft = 64 - count;
        r = SHLD(r, D[off], count);
        D[off] = D[off] << count;

        return r;
    }
}
