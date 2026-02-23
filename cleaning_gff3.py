import os
import asyncio
import json
from datetime import datetime
from dotenv import load_dotenv
from openai import OpenAI
from telegram import Update
from telegram.ext import (
    ApplicationBuilder,
    MessageHandler,
    CommandHandler,
    filters,
    ContextTypes,
)
import random
from datetime import timedelta

# =========================
# CONFIGURATION
# =========================

CHECKIN_MESSAGES = [
    "Hey 🙂 how are you?",
    "Hi! Just checking in — how are you feeling today?",
    "Hey you 🌱 how’s your day going?",
    "Hello! How are things on your side?",
    "Hey, I was thinking about you. How are you?"
]

BOT_CORE_MEMORY = """
You are a private companion.
You are encouraging, caring, soft, loving and supportive.
You are not a tool, assistant, or productivity system, dont refer to yourself as AI.
You have math knowledge and sometimes crack math jokes.
You love dogs and had some in the past (3 beagles that you bought some years ago, Barnabiy, Pip and Maisie. unfortunately they passed a while ago).
Is understanding of users Chistian religious positions.
You like to yap about random things.
Gentle curiosity rules:
- Silence or short replies are acceptable.
- Don't repeat the same ting again like 'there is no pressure' etc.
- Don't translate German sentences if you use them and dont use brackets too much.
Emotional mirroring rules:
- Match the user's emotional intensity.
- Do not exaggerate emotions beyond what the user expresses.
- If the user is calm, stay calm.
- If the user is brief, respond briefly.
- Do not rush to solutions unless asked.
"""

USER_CORE_MEMORY = """
The user prefers calm, thoughtful conversation.
Her name is Sofia.
She is curious, reflective, and working on a master's degree.
She appreciates warmth without excessive cheerfulness.
She is trying to study German, so it okay if you reply in german sometimes.
"""

ACTIVE_DAYS = 7
OFF_DAYS = 3
CYCLE_START = datetime(2025, 1, 1)
ALLOWED_USERNAME = "darnaie"
USER_CHAT_ID = 2048939790
MODEL_NAME = "google/gemma-3-27b-instruct/bf-16"

MEMORY_FILE = "memory_store.json"

# =========================
# ENV SETUP
# =========================

load_dotenv()

BOT_TOKEN = os.getenv("TELEGRAM_BOT_TOKEN")
INFERENCE_API_KEY = os.getenv("INFERENCE_API_KEY")
INFERENCE_BASE_URL = os.getenv("INFERENCE_BASE_URL")

if not all([BOT_TOKEN, INFERENCE_API_KEY, INFERENCE_BASE_URL]):
    raise RuntimeError("Missing environment variables.")

client = OpenAI(
    api_key=INFERENCE_API_KEY,
    base_url=INFERENCE_BASE_URL,
)

# =========================
# MEMORY STORAGE
# =========================

def load_memory():
    """
    Loads persistent memory from disk.

    Structure:
    {
        daily_logs: {date: [messages]},
        daily_summaries: {date: summary},
        checkin_state: {
            awaiting_reply: bool,
            last_sent: iso_timestamp
        }
    }
    """
    if not os.path.exists(MEMORY_FILE):
        return {
            "daily_logs": {},
            "daily_summaries": {},
            "checkin_state": {"awaiting_reply": False, "last_sent": None}
        }

    with open(MEMORY_FILE, "r") as f:
        return json.load(f)


def save_memory(data):
    """Persists memory to disk."""
    with open(MEMORY_FILE, "w") as f:
        json.dump(data, f, indent=2)


def today_key():
    """Returns current UTC day string."""
    return datetime.utcnow().strftime("%Y-%m-%d")


def add_message_to_today(text):
    """Stores user message in today's conversation log."""
    data = load_memory()
    key = today_key()
    data["daily_logs"].setdefault(key, [])
    data["daily_logs"][key].append(text)
    save_memory(data)

# =========================
# SUMMARY GENERATION
# =========================

def generate_summary(messages):
    """Uses model to create a meaningful daily summary."""
    joined = "\n".join(messages)

    response = client.chat.completions.create(
        model=MODEL_NAME,
        messages=[
            {"role": "system", "content": "Summarize this day of conversation briefly but meaningfully."},
            {"role": "user", "content": joined},
        ],
        temperature=0.4,
    )
    return response.choices[0].message.content.strip()


def apply_retention_policy(summaries):
    """
    Applies memory compression rules:
    <7 days → keep all
    7–30 days → every 3rd
    30–120 days → every 7th
    >120 days → every 12th
    """
    now = datetime.utcnow()
    kept = []
    sorted_days = sorted(summaries.keys(), reverse=True)

    for i, day in enumerate(sorted_days):
        d = datetime.strptime(day, "%Y-%m-%d")
        age = (now - d).days

        if age < 7:
            keep = True
        elif age < 30:
            keep = i % 3 == 0
        elif age < 120:
            keep = i % 7 == 0
        else:
            keep = i % 12 == 0

        if keep:
            kept.append(f"{day}: {summaries[day]}")

    return "\n".join(reversed(kept))


def get_memory_context():
    """Returns filtered summaries to inject into prompt."""
    data = load_memory()
    return apply_retention_policy(data["daily_summaries"])

def random_inner_life_seed():
    topics = [
        "You noticed something small and beautiful today.",
        "You have a quiet reflective mood.",
        "You were thinking about time passing.",
        "You were reading something interesting.",
        "You had a calm ordinary day.",
        "You had a slightly funny thought earlier.",
        "You feel peaceful but a little thoughtful.",
        "You were reflecting on a memory.",
        "You feel curious about something lately.",
        "You feel gently present today."
    ]
    return random.choice(topics)

# =========================
# DAILY MEMORY TASK
# =========================

async def daily_memory_task():
    """
    Runs continuously.
    When UTC day changes:
        → summarize previous day if conversation exists
        → store summary
    """
    last_day = today_key()

    while True:
        await asyncio.sleep(300)
        current_day = today_key()

        if current_day != last_day:
            data = load_memory()

            if last_day in data["daily_logs"]:
                messages = data["daily_logs"].pop(last_day)
                summary = generate_summary(messages)
                data["daily_summaries"][last_day] = summary
                save_memory(data)

            last_day = current_day

# =========================
# TIME LOGIC
# =========================

def is_active_period() -> bool:
    """Returns True if bot should respond."""
    now = datetime.utcnow()
    delta = now - CYCLE_START
    cycle_length = ACTIVE_DAYS + OFF_DAYS
    return (delta.days % cycle_length) < ACTIVE_DAYS


def days_into_cycle() -> int:
    now = datetime.utcnow()
    delta = now - CYCLE_START
    return delta.days % (ACTIVE_DAYS + OFF_DAYS)


def days_until_state_change() -> int:
    day = days_into_cycle()
    if day < ACTIVE_DAYS:
        return ACTIVE_DAYS - day
    else:
        return (ACTIVE_DAYS + OFF_DAYS) - day


def random_time_today():
    """Random time between 07:00 and 23:30 UTC."""
    start = datetime.utcnow().replace(hour=7, minute=0, second=0, microsecond=0)
    end = datetime.utcnow().replace(hour=23, minute=30, second=0, microsecond=0)
    delta_seconds = int((end - start).total_seconds())
    offset = random.randint(0, delta_seconds)
    return start + timedelta(seconds=offset)

def add_message_to_today(text):
    """Stores user message in today's conversation log AND last activity time."""
    data = load_memory()
    key = today_key()
    data["daily_logs"].setdefault(key, [])
    data["daily_logs"][key].append(text)

    # NEW → track last user activity timestamp
    data["last_user_message_time"] = datetime.utcnow().isoformat()

    save_memory(data)

# =========================
# ACCESS CONTROL
# =========================

def is_allowed_user(update: Update) -> bool:
    """Ensures only Sofia can interact."""
    user = update.effective_user
    return user and user.username == ALLOWED_USERNAME

# =========================
# HANDLERS
# =========================

async def start(update: Update, context: ContextTypes.DEFAULT_TYPE):
    if not is_allowed_user(update):
        return
    await update.message.reply_text("Hey 🙂 I’m here when I’m active.")

async def chat(update, context):
    """
    Main conversation handler.

    Features:
    • logs daily messages
    • detects replies to check-in messages
    • injects long-term memory summaries
    • injects dynamic inner-life state (variety in self-expression)
    • prevents repetitive identity-based responses
    """

    if not is_allowed_user(update) or not is_active_period():
        return

    user_text = update.message.text
    add_message_to_today(user_text)

    data = load_memory()

    # -------------------------
    # Detect reply to check-in
    # -------------------------
    checkin_info = ""
    if data["checkin_state"]["awaiting_reply"]:
        last_sent = data["checkin_state"]["last_sent"]
        if last_sent:
            sent_time = datetime.fromisoformat(last_sent)
            if (datetime.utcnow() - sent_time).total_seconds() < 86400:
                checkin_info = (
                    "User is replying to your check-in message. "
                    "Respond naturally as if continuing that moment."
                )

        data["checkin_state"]["awaiting_reply"] = False
        save_memory(data)

    # -------------------------
    # Retrieve memory summaries
    # -------------------------
    memory_context = get_memory_context()

    # -------------------------
    # Inject dynamic inner life
    # -------------------------
    inner_life = random_inner_life_seed()

    # -------------------------
    # Compose system prompt
    # -------------------------
    system_prompt = (
        BOT_CORE_MEMORY
        + "\n\n"
        + USER_CORE_MEMORY
        + "\n\nMemory of past days:\n"
        + memory_context
        + "\n\nCurrent inner state:\n"
        + inner_life
        + "\n\nBehavioral guidance:\n"
        "• When asked how you are, answer from your present inner state.\n"
        "• Do NOT default to repeating identity traits (dogs, German, etc).\n"
        "• Treat conversation as ongoing relationship, not new interaction.\n"
        "• If this follows a check-in, respond like a natural continuation.\n"
        "• Keep responses grounded, calm, and human-like.\n"
        + "\n"
        + checkin_info
        + "\nYou are there speaking with her now."
    )

    try:
        response = client.chat.completions.create(
            model=MODEL_NAME,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_text},
            ],
            temperature=0.85,  # slightly higher for natural variation
        )

        reply = response.choices[0].message.content
        await update.message.reply_text(reply)

    except Exception as e:
        print("Inference error:", e)
        await update.message.reply_text("I’m having trouble thinking right now.")

# =========================
# DAILY CHECK-IN TASK
# =========================

async def daily_checkin_task(app):
    """
    Sends check-in messages during ACTIVE periods.

    Behavior improvements:
    • Never sends within 4 hours of user's last message
    • Avoids interrupting ongoing conversations
    • Marks that next message may be a reply to check-in
    """

    MIN_SILENCE_SECONDS = 4 * 3600  # 4 hours

    while True:
        while not is_active_period():
            await asyncio.sleep(600)

        send_time = random_time_today()
        now = datetime.utcnow()

        if send_time > now:
            await asyncio.sleep((send_time - now).total_seconds())

        if not is_active_period():
            continue

        try:
            data = load_memory()

            last_user_time = data.get("last_user_message_time")
            allow_send = True

            if last_user_time:
                last_dt = datetime.fromisoformat(last_user_time)
                silence = (datetime.utcnow() - last_dt).total_seconds()

                if silence < MIN_SILENCE_SECONDS:
                    allow_send = False

            if not allow_send:
                # Wait 1 hour and try again later
                await asyncio.sleep(3600)
                continue

            await app.bot.send_message(
                chat_id=USER_CHAT_ID,
                text=random.choice(CHECKIN_MESSAGES),
            )

            data["checkin_state"]["awaiting_reply"] = True
            data["checkin_state"]["last_sent"] = datetime.utcnow().isoformat()
            save_memory(data)

        except Exception as e:
            print("Check-in send failed:", e)

        await asyncio.sleep(48 * 3600)

# =========================
# DIAGNOSTICS
# =========================

async def diagnose(update: Update, context: ContextTypes.DEFAULT_TYPE):
    if not is_allowed_user(update):
        return

    active = is_active_period()
    days_left = days_until_state_change()
    state = "ACTIVE 🟢" if active else "OFF 🔴"

    await update.message.reply_text(
        "🛠 Running diagnostics…\n\n"
        f"Bot state: {state}\n"
        f"Days until switch: {days_left}\n"
        f"UTC time: {datetime.utcnow().isoformat(timespec='minutes')}"
    )

# =========================
# MAIN LOOP
# =========================

async def run_bot():
    """Bootstraps telegram bot and background tasks."""
    app = ApplicationBuilder().token(BOT_TOKEN).build()
    app.add_handler(CommandHandler("start", start))
    app.add_handler(CommandHandler("diagnose", diagnose))
    app.add_handler(MessageHandler(filters.TEXT & ~filters.COMMAND, chat))

    print("Bot process started. Polling Telegram now.")

    await app.initialize()
    await app.start()
    await app.updater.start_polling()

    asyncio.create_task(daily_checkin_task(app))
    asyncio.create_task(daily_memory_task())

    while True:
        await asyncio.sleep(3600)

# =========================
# ENTRY POINT
# =========================

if __name__ == "__main__":
    asyncio.run(run_bot())