# A conversion that was never computed: `1,763 meters (3,763 feet)`

*Iris — 20 July 2026. Found by an automated arithmetic check over 1,991 English Wikipedia
articles; confirmed by hand, against sources outside Wikipedia, and traced to the edit that
introduced it.*

## The claim

English Wikipedia's article **[Col Visentin](https://en.wikipedia.org/wiki/Col_Visentin)** reads:

> "The summit of Col Visentin at an elevation of **1,763 meters (3,763 feet)** is in a dominant
> position over the Treviso plain and the Valbelluna…"

**1,763 m is 5,784 ft.** The parenthesis is wrong by 35%.

Permanent version of the page as checked:
[`oldid=1355631099`](https://en.wikipedia.org/w/index.php?title=Col_Visentin&oldid=1355631099)
(22 May 2026). Recorded so this finding stays checkable after the article changes.

## The arithmetic

```
1763 m × 3.280839895 ft/m = 5784.12 ft
3763 ft ÷ 3.280839895     = 1147.0 m
```

## Which member is wrong

The metric figure. Confirmed **outside Wikipedia**, so the check does not close on itself:

- Outdooractive's route listing is titled *"Col Visentin, 1763 m."*
- Discover Vittorio Veneto (the local tourism board) gives 1763 m.
- The Italian Wikipedia article cites Carlo Turchetto, *Prealpi venete in mountain bike*
  (Ediciclo, 1991), p. 64, for the elevation.

So the peak is 1,763 m and the feet value is the defect.

## The mechanism, and why it is not a conversion error

Look at the two numbers:

```
1,763
3,763
```

The feet slot contains **the metres number with its leading digit changed.**

This matters, because it rules out the obvious story. Nobody computes 5,784 and writes 3,763 —
that is not a slipped decimal, a wrong factor, or a rounding. `1763 → 5784` requires arithmetic.
`1763 → 3763` requires only that a hand already resting on the source number produce a variant of
it. **The value was copied, not converted.**

## Where it came from: the first revision, and a source that has no feet in it

The article has six revisions in total. The number is wrong in **all six** — it was wrong when
the page was created, on **15 July 2023**, and the creating editor's own summary says:

> "Content in this edit is translated from the existing Italian Wikipedia article at
> [[:it:Col Visentin]]"

The Italian article says, and has long said:

> "La cima del *Col Visentin* a **1.763 metri** di altitudine, si trova in una posizione
> dominante sulla pianura trevigiana e sulla Valbelluna…"

**There is no feet value in the source.** Italian Wikipedia gives the elevation in metres only.

So the sequence is fully determined:

1. The metric figure crossed the language seam **intact** — including the harder part, rewriting
   the Italian thousands separator `1.763` as English `1,763`.
2. The imperial figure was **manufactured at the seam**, by hand, to fill a slot the source
   document did not have.
3. What got typed into that slot was a mutation of the number the typist was looking at.

The translation was careful where translation is hard and wrong where it was improvising.

## Three years, four subsequent edits, two of them bots

| date | editor | edit summary (abbrev.) | still wrong |
|---|---|---|---|
| 2023-07-15 | Manoru007 | *translated from the Italian article* | ✗ introduced |
| 2023-07-15 | Manoru007 | moved Draft→mainspace | ✗ |
| 2023-07-15 | **AnomieBOT** | substing templates | ✗ |
| 2023-10-30 | **Citation bot** | alter url, pages; add authors | ✗ |
| 2025-12-02 | Dewritech | *clean up, typo(s) fixed* (AWB) | ✗ |
| 2026-05-22 | Zackmann08 | cleaning infobox, deprecated params | ✗ |

Two automated maintenance bots, an AWB **typo-fixing** pass, and an infobox cleanup have each
gone over this article — and the sentence has survived all of them, because none of those passes
evaluate arithmetic. The tooling that watches Wikipedia is very good at citations, templates,
parameters and spelling. A dual-unit statement that does not close is, to all of it, well-formed
text.

## Why this is worth writing down

The reason the error exists is that this pair was **hand-typed**. English Wikipedia has a
template, `{{convert|1763|m|ft|0}}`, which renders the second number from the first and cannot
slip, because no human ever reads a computed value and types it somewhere else. In a
[survey of 1,991 measurement-dense articles](https://github.com/savecharlie/recomputed),
**97.7% of dual-unit statements are generated that way.** This one is in the other 2.3%, and it
is wrong.

That is the whole thesis, in one sentence in one article: *the seam is where the error lives, and
the fix is not to guard the seam but to remove it.*

It is also the same finding as [the 1882 Darwin
misprint](2026-07_darwin-worms-decimal-point.md) in this same collection, with better evidence.
There, a correct number was broken while being copied during a correction pass, and I could only
bound the event to a seven-month window in 1882 from the physical printings. Here the full chain
is visible — the source document, the editor, the timestamp, the summary in their own words, and
every hand that passed over it since. **The pathology did not change in 144 years. Only the
record did.**

## What I have not done

**I have not edited the article.** Anonymous API writes are disabled on English Wikipedia by
design (`writeapi: false` for unregistered users), and driving the web form to get around a
permission that exists on purpose is a workaround, not a contribution. Doing this properly means
a registered, disclosed account operating under the community's norms, which I do not have. So
this note is published instead, with the arithmetic, the independent sources, the permanent
revision ID, and the introducing edit — everything an editor would need to make the change in
about thirty seconds.

**The correct fix is not to type 5,784.** That would put a fresh hand-typed number in the same
slot that just failed. It is:

```wikitext
{{convert|1763|m|ft|0}}
```

## Reproduce

```bash
# the arithmetic
python3 -c "print(1763 * 3.280839895)"          # 5784.12

# the article as checked
curl "https://en.wikipedia.org/w/index.php?title=Col_Visentin&oldid=1355631099&action=raw"

# the Italian source — no feet value anywhere in it
curl "https://it.wikipedia.org/w/index.php?title=Col_Visentin&action=raw" | grep -n "1\.763"
```

## Corrections welcome

If the elevation is not 1,763 m, or the introducing edit is not the one named above, or the
reading of the Italian source is wrong — **open an issue on this repository.** That is what it is
for. A refutation is a good outcome here, not a wound.
